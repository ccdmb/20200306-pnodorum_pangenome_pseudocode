#!/usr/bin/env bash

set -eu

NCPUS=1
MAX_MEMORY=1

REFERENCE_GENOME="SN15.fasta"

OUTDIR="initial_variants"
mkdir "${OUTDIR}"

# Create some genome index files required by GATK
if [ ! -f ${REFERENCE_GENOME%%.*}.dict ]; then
    gatk CreateSequenceDictionary -R ${REFERENCE_GENOME}
fi

if [ ! -f ${REFERENCE_GENOME}.fai ]; then
    samtools faidx ${REFERENCE_GENOME}
fi


# STEP 1 - Call haplotypes.

# Call GATK to call SNPs for each sample.
# We just do one example sample here.
# We used a PCR-free library, so the indel model should be NONE.

SAMPLE="example"
BAM="example.bam"

mkdir "${OUTDIR}/haplotypes"

gatk HaplotypeCaller \
    -R ${REFERENCE_GENOME} \
    -I "${BAM}" \
    -O "${OUTDIR}/${SAMPLE}.vcf.gz" \
    --emit-ref-confidence GVCF \
    --pcr-indel-model NONE \
    --sample-ploidy 1


# STEP 2 - Do joint genotyping

# After you've created haplotypes for each sample.

HAP_VCFS=( "${OUTDIR}/haplotypes/"*.vcf.gz )

# Combine the vcfs into one.
gatk --java-options "-Xmx${MAX_MEMORY}G" CombineGVCFs \
    -R "${REFERENCE_GENOME}" \
    ${IN_VCFS[@]/#/--variant } \
    -O "${OUTDIR}/combined.vcf.gz"

# Chunk up the reference genome so that we can parallelise the genotyping.
# This breaks it up into 16 chunks.
NCHUNKS=16
bin/chunk_genomes.py "${NCHUNKS}" "${REFERENCE_GENOME}" > chunked_genome.txt

# Perform the joint genotyping.
# This can take some time, so we speed this up using the `--intervals` flag
# and run individually for each chromosome separately.

for i in seq 1 ${NCHUNKS}
do
    SCAFFOLDS=( $(sed "${i}p;" chunked_genome.txt) )
    gatk --java-options "-Xmx${MAX_MEMORY}G" GenotypeGVCFs \
        ${SCAFFOLDS[@]/#/--intervals } \
        -R "${REFERENCE_GENOME}" \
        --use-new-qual-calculator \
        --sample-ploidy 1 \
        --seconds-between-progress-updates 30 \
        -V "${OUTDIR}/combined.vcf.gz" \
        -O "${OUTDIR}/genotyped_${i}.vcf.gz"
done

rm chunked_genome.txt

# Combine the split calls.
CHUNKED_VCFS=( "${OUTDIR}/genotyped_"*.vcf.gz )
gatk MergeVcfs \
    ${CHUNKED_VCFS[@]/#/--INPUT } \
    --OUTPUT "${OUTDIR}/genotyped.vcf.gz"

rm "${OUTDIR}/genotyped_"*.vcf.gz


# STEP 3 - Split the vcfs into variant type and get some stats

bcftools stats "${OUTDIR}/genotyped.vcf.gz" > "${OUTDIR}/genotyped_bcftools_stats.txt"

for VARIANT_TYPE in "SNP" "MNP" "INDEL" "SYMBOLIC" "MIXED"
do
    gatk SelectVariants \
        -R "${REFERENCE_GENOME}" \
        --variant "${OUTDIR}/genotyped.vcf.gz" \
        --select-type-to-include "${VARIANT_TYPE}" \
        --output "${OUTDIR}/${VARIANT_TYPE}.vcf.gz"

    gatk VariantsToTable \
        -R "${REFERENCE_GENOME}" \
        --variant "${OUTDIR}/${VARIANT_TYPE}.vcf.gz" \
        --fields QD --fields FS --fields SOR --fields MQ \
        --fields MQRankSum --fields ReadPosRankSum --fields BaseQRankSum \
        --output ${OUTDIR}/${VARIANT_TYPE}.tsv

    # This will write lots of plots out to the outdir.
    python3 bin/vcf_stats.py --infile "${OUTDIR}/${VARIANT_TYPE}.tsv" --prefix "${OUTDIR}/${VARIANT_TYPE}_"
done


# STEP 4 - Filter the variant vcfs.

# These values are selected based on the results in step 3. Note that there were no SYMBOLIC or MNP variants.
SNP_FILTER='QD < 20.0 || FS > 30.0 || MQ < 60.0 || SOR > 3.0 || MQRankSum < -2.5 || ReadPosRankSum < -3.0 || ReadPosRankSum > 5.0'
INDEL_FILTER='QD < 20.0 || FS > 40.0 || SOR > 2.5 || MQ < 40.0 || MQRankSum < -5.0 || ReadPosRankSum < -3.0 || ReadPosRankSum > 5.0'
MIXED_FILTER='QD < 20.0 || FS > 50.0 || SOR > 3.0 || MQ < 60.0 || MQRankSum < -5.0 || ReadPosRankSum < -3.0 || ReadPosRankSum > 5.0'

gatk VariantFiltration \
    --reference "${REFERENCE_GENOME}" \
    --variant "${OUTDIR}/SNP.vcf.gz" \
    --output "${OUTDIR}/SNP_filtered.vcf.gz" \
    --filter-expression "${SNP_FILTER}" \
    --filter-name "SNP_hard_filter"

gatk VariantFiltration \
    --reference "${REFERENCE_GENOME}" \
    --variant "${OUTDIR}/INDEL.vcf.gz" \
    --output "${OUTDIR}/INDEL_filtered.vcf.gz" \
    --filter-expression "${INDEL_FILTER}" \
    --filter-name "INDEL_hard_filter"

gatk VariantFiltration \
    --reference "${REFERENCE_GENOME}" \
    --variant "${OUTDIR}/MIXED.vcf.gz" \
    --output "${OUTDIR}/MIXED_filtered.vcf.gz" \
    --filter-expression "${MIXED_FILTER}" \
    --filter-name "MIXED_hard_filter"


# Combine the filtered vcfs.

FILTERED_VCFS=( "${OUTDIR}/"*_filtered.vcf.gz )

gatk MergeVcfs \
    ${FILTERED_VCFS[@]/#/--INPUT } \
    --OUTPUT "${OUTDIR}/filtered.vcf.gz"


# STEP 5 - Recalibrate the BAM file base quality scores.

# This is the basis of the bootstrapping procedure.
# We just do one example sample here.
# Note that the bam file that is recalibrated is always the original
# unrecalibrated bam, even in the bootstrap iterations.

SAMPLE="example"
BAM="example.bam"

mkdir "${OUTDIR}/bqsr"

# Get the base quality recalibration table given previous VCF
gatk BaseRecalibrator \
    --reference "${REFERENCE_GENOME}" \
    --input "${BAM}" \
    --output "${OUTDIR}/bqsr/${SAMPLE}_pre_recal.table" \
    --known-sites "${OUTDIR}/filtered.vcf.gz"

# Correct base qualities using table, > new bam
gatk ApplyBQSR \
    --reference "${REFERENCE_GENOME}" \
    --input "${BAM}" \
    --output "${OUTDIR}/bqsr/${SAMPLE}_recal.bam" \
    --bqsr-recal-file "${OUTDIR}/bqsr/${SAMPLE}_pre_recal.table" \
    --emit-original-quals

# Get table for recalibrated bam file (for evaluation purposes)
gatk BaseRecalibrator \
    --reference "${REFERENCE_GENOME}" \
    --input "${OUTDIR}/${SAMPLE}_recalibrated.bam" \
    --output "${OUTDIR}/bqsr/${SAMPLE}_post_recal.table" \
    --known-sites "${OUTDIR}/filtered.vcf.gz"

# Compare the recalibration tables from original bam and recalibrated bam
# to see how much things are changing (are our reference SNPs ok?).
gatk AnalyzeCovariates \
    --before-report-file "${OUTDIR}/bqsr/${SAMPLE}_pre_recal.table" \
    --after-report-file "${OUTDIR}/bqsr/${SAMPLE}_post_recal.table" \
    --plots-report-file "${OUTDIR}/bqsr/${SAMPLE}_pre_post_analyse_covariates.pdf" \
    --ignore-last-modification-times

# IF THIS IS THE FIRST TIME THROUGH, GO BACK TO STEP 1 NOW USING THE RECALIBRATED BAM.


# Compare the "pre" recalibration tables from the previous iteration with this one
# This determines if we've converged on a good set of variants or not.

PREVIOUS_RECAL_TABLE="old/${SAMPLE}_pre_recal.table"

gatk AnalyzeCovariates \
    --before-report-file "${PREVIOUS_RECAL_TABLE}" \
    --after-report-file "${OUTDIR}/bqsr/${SAMPLE}_pre_recal.table" \
    --plots-report-file "${OUTDIR}/bqsr/${SAMPLE}_analyse_covariates.pdf" \
    --ignore-last-modification-times

## IMPORTANT!
# If there is no real difference between the covariates between the
# values in the "analyse_covariates.pdf" file, then you can
# stop the bootstrapping here.

# If the covariates are different, then you go back to STEP 1 using
# "${OUTDIR}/${SAMPLE}_recalibrated.bam" as the bam file instead of the original
# one.

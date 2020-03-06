#!/usr/bin/env bash

set -eu

NCPUS=1
MAX_MEMORY=1

REFERENCE_GENOME="SN15.fasta"
REFERENCE_GENES="SN15.gff3"

OUTDIR="filtered_variants"
VCF="combined.vcf.gz"

# Get a .fai file.
samtools faidx "${REFERENCE_GENOME}"

mkdir "${OUTDIR}"


# STEP 1 - Get statistics for individual sample genotypes in the vcf file.

CHROMS=$(cut -f1 "${REFERENCE_GENOME}.fai")

for CHROM in ${CHROMS}
do
    bcftools view \
      --regions "${CHROM}" \
      --threads "${NCPUS}" \
      -l 1 \
      -O z \
      -o "${OUTDIR}/${CHROM}.vcf.gz" \
      "${VCF}"

    bin/indiv_snps.R \
      --vcf "${OUTDIR}/${CHROM}.vcf.gz" \
      --outdir "${OUTDIR}/${CHROM}" \
      --no-prefilter \
      --no-postfilter \
      --min-gq 0.3 \
      --min-dp 8 \
      2>&1 \
    | tee "${OUTDIR}/${CHROM}.log"

    rm "${OUTDIR}/${CHROM}.vcf.gz"
done


# STEP 2 - Apply genotype filters determined in step 1.

MIN_DP=8
MIN_GQ=30

gatk VariantFiltration \
  --genotype-filter-expression "DP < ${MIN_DP}" \
  --genotype-filter-name "minGTDP" \
  --variant "${VCF}" \
  --output "${OUTDIR}/gt_filtered.tmp1.vcf.gz" \
  --set-filtered-genotype-to-no-call \
  --reference "${REFERENCE_GENOME}"

gatk VariantFiltration \
  --genotype-filter-expression "GQ < ${MIN_GQ}" \
  --genotype-filter-name "minGTGQ" \
  --variant "${OUTDIR}/gt_filtered.tmp1.vcf.gz" \
  --output "${OUTDIR}/gt_filtered.tmp2.vcf.gz" \
  --set-filtered-genotype-to-no-call \
  --reference "${REFERENCE_GENOME}"

bcftools view \
  --trim-alt-alleles \
  -Ou \
  ${OUTDIR}/gt_filtered.tmp2.vcf.gz \
| bcftools filter \
  --exclude "F_PASS(GT='mis' | GT='ref') == 1" \
  --soft-filter AllGTRefOrMissing \
  -O z \
  -o "${OUTDIR}/gt_filtered.tmp3.vcf.gz"

bcftools index -t "${OUTDIR}/gt_filtered.tmp3.vcf.gz"

gatk LeftAlignAndTrimVariants \
  --reference "${REFERENCE_GENOME}" \
  -V "${OUTDIR}/gt_filtered.tmp3.vcf.gz" \
  -O "${OUTDIR}/gt_filtered.vcf.gz"


rm "${OUTDIR}"/gt_filtered.tmp*.vcf.gz "${OUTDIR}"/gt_filtered.tmp*.vcf.gz.tbi


# STEP 3 - Visualise variant loci from the genotype filtered set.

VARIANT_TYPES="SNP MNP INDEL SYMBOLIC MIXED"

for TYPE in ${VARIANT_TYPES}
do

    mkdir -p "${OUTDIR}/${TYPE}"

    gatk SelectVariants \
      -R "${GENOME}" \
      --variant "${INFILE}" \
      --select-type-to-include ${TYPE} \
      --output "${OUTDIR}/${TYPE}.vcf.gz"

    gatk VariantsToTable \
      -R ${REFERENCE_GENOME} \
      --variant "${OUTDIR}/${TYPE}.vcf.gz" \
      --fields QD --fields FS --fields SOR --fields MQ \
      --fields MQRankSum --fields ReadPosRankSum --fields BaseQRankSum \
      --output "${OUTDIR}/${TYPE}.tsv"

    python3 bin/vcf_stats.py \
      --infile "${OUTDIR}/${TYPE}.tsv" \
      --prefix "${OUTDIR}/${tYPE}/"

    bin/chrom_snps.R \
      --vcf "${OUTDIR}/${TYPE}.vcf.gz" \
      --fasta "${REFERENCE_GENOME}" \
      --gff "${REFERENCE_GENES}" \
      --outdir "${OUTDIR}/${TYPE}" \
      --max-dp 12000
done


# STEP 4 - Filter genotype loci using thresholds determined in step 3.

read -d '' TABLE <<EOF || true
SNP,QD,QD < 15.0
SNP,FS,FS > 60.0
SNP,SOR,SOR > 3.5
SNP,MQ,MQ < 70.0
SNP,MQRanksum,MQRankSum < -5.0
SNP,ReadPosRankSum,ReadPosRankSum < -2.5 || ReadPosRankSum > 3.0
SNP,DP,DP < 2000 || DP > 9500
INDEL,QD,QD < 15.0
INDEL,FS,FS > 50.0
INDEL,SOR,SOR > 4.0
INDEL,MQ,MQ < 70.0
INDEL,MQRankSum,MQRankSum < -6.0
INDEL,ReadPosRankSum,ReadPosRankSum < -4.0 || ReadPosRankSum > 3.5
INDEL,DP,DP < 2500 || DP > 9000
MIXED,QD,QD < 15.0
MIXED,FS,FS > 50.0
MIXED,SOR,SOR > 3.5
MIXED,MQ,MQ < 70.0
MIXED,MQRankSum,MQRankSum < -5.0
MIXED,ReadPosRankSum,ReadPosRankSum < -3.0 || ReadPosRankSum > 3.0
MIXED,DP,DP < 2500 || DP > 9000
MNP,QD,QD < 15.0
MNP,FS,FS > 50.0
MNP,SOR,SOR > 3.5
MNP,MQ,MQ < 70.0
MNP,MQRankSum,MQRankSum < -5.0
MNP,ReadPosRankSum,ReadPosRankSum < -3.0 || ReadPosRankSum > 3.0
MNP,DP,DP < 2500 || DP > 9000
EOF

while read LINE; do
  TYPE="$(echo ${LINE} | cut -d, -f 1)"
  STAT="$(echo ${LINE} | cut -d, -f 2)"
  FILTER="$(echo ${LINE} | cut -d, -f 3)"

  gatk VariantFiltration \
    --reference ${GENOME} \
    --variant "${OUTDIR}/${TYPE}.vcf.gz" \
    --output "${OUTDIR}/${TYPE}.tmp.vcf.gz" \
    --filter-expression "${FILTER}" \
    --filter-name "${TYPE}_${STAT}"

  mv "${OUTDIR}/${TYPE}.tmp.vcf.gz" "${OUTDIR}/${TYPE}.vcf.gz"
  mv "${OUTDIR}/${TYPE}.tmp.vcf.gz.tbi" "${OUTDIR}/${TYPE}.vcf.gz.tbi"

done < <(echo -e "${TABLE}")

VTYPES=( $(echo -e "${TABLE}" | cut -d, -f 1 | uniq) )
VTYPES=( ${VTYPES[@]/%/.vcf.gz})
VTYPES=( ${VTYPES[@]/#/${OUTDIR}/})

gatk MergeVcfs \
  ${VTYPES[@]/#/-I } \
  -R ${GENOME} \
  -O ${OUTDIR}/combined.vcf.gz


# Run steps 1 and 3 again using the filtered variants.

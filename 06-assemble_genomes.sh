#!/usr/bin/env bash

set -eu

NCPUS=1
MAX_MEMORY=4

GENOME=SN15.fasta
GENES=SN15.gff3

FASTQ_DIR=fastq
OUTDIR=06-assemble_genomes

BASE_NAME="example"
R1_FQ=( "example_pair1_R1.fastq.gz" "example_pair2_R1.fastq.gz" )
R2_FQ=( "example_pair1_R2.fastq.gz" "example_pair2_R2.fastq.gz" )
SAMPLES=( "example_pair1" "example_pair2" )


# STEP 1 - Stitch paired end reads for samples with shorter insert size.
mkdir -p "${OUTDIR}/fastq"

bbmerge-auto.sh \
  in1="fastq/${R1_FQ[0]}" \
  in2="fastq/${R2_FQ[0]}" \
  out="${OUTDIR}/fastq/${SAMPLES[0]}" \
  outu1="${OUTDIR}/fastq/${R1_FQ[0]}" \
  outu2="${OUTDIR}/fastq/${R2_FQ[0]}" \
  ihist="${OUTDIR}/fastq/${SAMPLES[0]}_insert.txt" \
  vstrict \
  rem \
  k=62 \
  extend2=50 \
  ecct


# STEP 2 - Assemble the genomes
mkdir -p "${OUTDIR}/assemblies"

# Formats the read files to the parameter format required for spades.
format_read_args () {
    ARR=( $1 )
    READ=$2
    STR=""
    for i in "${!ARR[@]}"; do
        j=$(( i + 1 ))
        STR+=" --pe${j}-${READ} ${ARR[i]}"
    done
    echo ${STR}
}

## STEP 2a - Run assembly for non-stitched samples.

# Actually calls the arg formatter.
FWD_ARGS=$(format_read_args "${R1_FQ[@]/#/${FASTQ_DIR}/}" 1)
REV_ARGS=$(format_read_args "${R2_FQ[@]/#/${FASTQ_DIR}/}" 2)

spades.py \
    -o "${OUTDIR}/assemblies/${BASE_NAME}" \
    ${FWD_ARGS} \
    ${REV_ARGS} \
    --careful -k 21,31,51,71,81,101 \
    --cov-cutoff auto \
    --memory ${MAX_MEMORY} \
    --threads ${NCPUS}

## STEP 2b - Run assembly for non-stitched samples.

# Actually calls the arg formatter.
FWD_ARGS=$(format_read_args "${R1_FQ[@]/#/${OUTDIR}/fastq/}" 1)
REV_ARGS=$(format_read_args "${R2_FQ[@]/#/${OUTDIR}/fastq/}" 2)
MERGED_ARGS=$(format_read_args "${R2_FQ[@]/#/${OUTDIR}/fastq/}" "m")

spades.py \
    -o "${OUTDIR}/${BASE_NAME}" \
    ${FWD_ARGS} \
    ${REV_ARGS} \
    --careful -k 31,51,71,81,101,127 \
    --cov-cutoff auto \
    --memory ${MAX_MEMORY} \
    --threads ${NCPUS}


# STEP 3 - Run quast on the assemblies
mkdir -p "${OUTDIR}/quast"

## STEP 3a - Run quast for non-stitched samples.
# Formats the read files to the parameter format required for spades.
format_read_args () {
    ARR1=( $1 )
    ARR2=( $2 )
    STR=""
    for i in "${!ARR1[@]}"; do
        j=$(( i + 1 ))
        STR+=" --pe1 ${ARR1[i]}"
        STR+=" --pe2 ${ARR2[i]}"
    done
    echo ${STR}
}

quast.py \
    -R "${GENOME}" \
    -g "${GENES}" \
    -o "${OUTDIR}/quast/${SAMPLE}" \
    --min-contig 500 \
    --rna-finding \
    --threads ${NCPUS} \
    --gene-finding \
    --conserved-genes-finding \
    --fungus \
    --labels "${SAMPLE}_scaffolds ${SAMPLE}_contigs" \
    $(format_read_args "${R1_FQ[@]/#/${FASTQ_DIR}/}" "${R2_FQ[@]/#/${FASTQ_DIR}/}") \
    "${OUTDIR}/assemblies/${SAMPLE}/scaffolds.fasta" \
    "${OUTDIR}/assemblies/${SAMPLE}/contigs.fasta"


## STEP 3b - Run quast for the stitched samples

# Formats the read files to the parameter format required for spades.
format_read_args () {
    ARR1=( $1 )
    ARR2=( $2 )
    ARRM=( $3 )
    STR=""
    for i in "${!ARR1[@]}"; do
        j=$(( i + 1 ))
        STR+=" --pe1 ${STEP_MERGE}/${ARR1[i]}"
        STR+=" --pe2 ${STEP_MERGE}/${ARR2[i]}"
        STR+=" --single ${STEP_MERGE}/${ARRM[i]}"
    done
    echo ${STR}
}

quast.py \
    -R "${GENOME}" \
    -g "${GENES}" \
    -o "${OUTDIR}/quast/${SAMPLE}" \
    --min-contig 500 \
    --rna-finding \
    --threads ${NCPUS} \
    --gene-finding \
    --conserved-genes-finding \
    --fungus \
    --labels "${SAMPLE}_scaffolds ${SAMPLE}_contigs" \
    $(format_read_args "${R1_FQ[@]/#/${OUTDIR}/fastq/}" "${R2_FQ[@]/#/${OUTDIR}/fastq/}" "${SAMPLES[@]/#/${OUTDIR}/fastq/}") \
    "${OUTDIR}/assemblies/${SAMPLE}/scaffolds.fasta" \
    "${OUTDIR}/assemblies/${SAMPLE}/contigs.fasta"

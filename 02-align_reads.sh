#!/usr/bin/env bash

set -eu

mkdir stampy_index

NCPUS=1
MAX_MEMORY=1

PICARD="/path/to/picard.jar"

# These are all dummy names.
REFERENCE_GENOME="SN15.fasta"

BASE_NAME="example"
R1_FQ=( "example_pair1_R1.fastq.gz" "example_pair2_R1.fastq.gz" )
R2_FQ=( "example_pair1_R2.fastq.gz" "example_pair2_R2.fastq.gz" )
SAMPLES=( "example" "example" )

# This is composed of the flowcell id, the flowcell lane number, and the sample name.
READ_GROUPS=( "HN2KYCCXY.7.example" "HTH77CCXY.8.example" )

PLATFORM="illumina"

# STEP 1 - Index the reference genome

stampy.py \
  --species="Pnodorum" \
  --assembly="SN15" \
  --build-genome=stampy_index/index \
  "${REFERENCE_GENOME}"

stampy.py \
    --genome=stampy_index/index \
    --build-hash=stampy_index/hash


# STEP 2 - Align the reads!
# You should do this for all read pairs.

PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"
REV_READ="${R2_FQ[${PAIR_INDEX}]}"
SAMPLE="${SAMPLE[${PAIR_INDEX}]}"
READ_GROUP="${READ_GROUPS[${PAIR_INDEX}]}"

stampy.py \
    -g stampy_index/index \
    -h stampy_index/hash \
    --threads "${NCPUS}" \
    --sensitive \
    --substitutionrate=0.0001 \
    --readgroup="ID:${READ_GROUP},PU:${READ_GROUP},SM:${SAMPLE},LB:${SAMPLE},PL:${PLATFORM}" \
    --xa-max=3 \
    --xa-max-discordant=10 \
    --outputformat=sam \
    --output="${SAMPLE}_pair${PAIR_INDEX}.sam" \
    -M "${FWD_READ}" "${REV_READ}"


samtools sort \
    -l 9 \
    -O bam \
    -@ "${NCPUS}" \
    "${SAMPLE}_pair${PAIR_INDEX}.sam" \
> "${SAMPLE}_pair${PAIR_INDEX}.bam"

# Remove the intermediate SAM file
rm "${SAMPLE}_pair${PAIR_INDEX}.sam"

# Calculate the bam index file.
samtools index "${SAMPLE}_pair${PAIR_INDEX}.bam"

# Calculate some basic stats
fastqc --threads "${NCPUS}" "${SAMPLE}_pair${PAIR_INDEX}.bam"

samtools stats \
    --remove-dups \
    -@ "${NCPUS}" \
    "${SAMPLE}_pair${PAIR_INDEX}.bam" \
> "${SAMPLE}_pair${PAIR_INDEX}.bam.stats.txt"

samtools flagstat \
    -@ "${NCPUS}" \
    "${SAMPLE}_pair${PAIR_INDEX}.bam" \
> "${SAMPLE}_pair${PAIR_INDEX}.bam.flagstat.txt"



# STEP 3 - Combine all of the fastq files and mark duplicates.

# After the last step, we'll have multiple bams like this:
ALIGNED_BAMS=( "${SAMPLE}_pair0.bam" "${SAMPLE}_pair1.bam" )

# The read groups are all marked in the SAM fields, so we can merge
# them to form a single bam per sample.

java -Xmx${MAX_MEMORY}G -jar "${PICARD}" MarkDuplicates \
    ${ALIGNED_BAMS[@]/#/INPUT=} \
    OUTPUT="${SAMPLE}.bam" \
    METRICS_FILE="${SAMPLE}_mark_duplicates.txt"

samtools index "${SAMPLE}.bam"

#!/usr/bin/env bash

set -eu

NCPUS=1
MAX_MEMORY=4

GENOME=SN15.fasta
MITO=SN15_mito.fasta
GENES=SN15.gff3

KMER=39
MIN_SIZE=12000
MAX_SIZE=120000
READ_LENGTH=125 # 150
INSERT_SIZE=600 # 300

BASE_NAME="example"
R1_FQ=( "example_pair1_R1.fastq.gz" "example_pair2_R1.fastq.gz" )
R2_FQ=( "example_pair1_R2.fastq.gz" "example_pair2_R2.fastq.gz" )
SAMPLES=( "example_pair1" "example_pair2" )

FASTQ_DIR=fastq
NUCLEAR_ASSEMBLY_DIR="06-assemble_genomes/assemblies/${BASENAME}"

OUTDIR=07-assemble_mitochondrial_genomes

mkdir -p "${OUTDIR}/mitoasm"

# STEP 1 - Assemble mito genomes

cat ${R1_FQ[@]/#/${FASTQ_DIR}/} > forward.fastq.gz
cat ${R2_FQ[@]/#/${FASTQ_DIR}/} > reverse.fastq.gz


cat << EOF > config.txt
Project:
-----------------------
Project name          = ${BASE_NAME}
Type                  = mito
Genome Range          = ${MIN_SIZE}-${MAX_SIZE}
K-mer                 = ${KMER}
Extended log          = 0
Save assembled reads  = no
Seed Input            = ${MITO}
Reference sequence    = ${MITO}
Variance detection    = no
Dataset 1:
-----------------------
Read Length           = ${READ_LENGTH}
Insert size           = ${INSERT_SIZE}
Platform              = illumina
Single/Paired         = PE
Forward reads         = forward.fastq.gz
Reverse reads         = reverse.fastq.gz
Optional:
-----------------------
Insert size auto      = yes
Insert Range          = 1.6
Insert Range strict   = 1.2
Use Quality Scores    = no
EOF

NOVOPlasty.pl -c config.txt

# Rename for better sorting and consistency.
if [[ -f  "Circularized_assembly_1_${NAME}.fasta" ]]; then
    mv "Circularized_assembly_1_${NAME}.fasta" "${OUTDIR}/mitoasm/${NAME}_mitochondrial.fasta"
elif [[ -f  Uncircularized_assemblies_1_${NAME}.fasta ]]; then
    mv "Uncircularized_assemblies_1_${NAME}.fasta" "${OUTDIR}/mitoasm/${NAME}_mitochondrial.fasta"
fi
mv "log_${NAME}.txt" "${OUTDIR}/mitoasm/${NAME}_log.txt"


# STEP 2 - Align scaffolds to mitochondria

mkdir -p "${OUTDIR}/alignments"

minimap2 \
    -x asm20 \
    -o "${OUTDIR}/alignments/${NAME}.paf \
    "${OUTDIR}/mitoasm/${NAME}_mitochondrial.fasta" \
    "${NUCLEAR_ASSEMBLY_DIR}/${NAME}/scaffolds.fasta"

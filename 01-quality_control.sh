#!/usr/bin/env bash

set -eu

MAX_MEMORY=1
NCPUS=1

TRIM_PHRED=2  # Do quality trimming on ends of sections below this score.
FILTER_PHRED=5  # Exclude reads with an average score below this.
MIN_READ_LENGTH=50  # The minimum length a read can be after trimming.


DATA_DIR="./data"  # Things like adapter sequences are stored in here.

# We'll store our results in these directories.
mkdir stats
mkdir trimmed
mkdir filtered_synthetic
mkdir kraken
mkdir alignments

# These are all dummy names.
BASE_NAME="example"
R1_FQ=( "example_pair1_R1.fastq.gz" "example_pair2_R1.fastq.gz" )
R2_FQ=( "example_pair1_R2.fastq.gz" "example_pair2_R2.fastq.gz" )


# STEP 1. Get raw statistics for each readpair.

## We're just doing one pair in the array as an example.
## This could be a for-loop or parallelised.
PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"
REV_READ="${R2_FQ[${PAIR_INDEX}]}"


bbcountunique.sh \
  -Xmx${MAX_MEMORY}g \
  in1="${FWD_READ}" \
  in2="${REV_READ}" \
  out="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_count_unique.txt"

kmercountmulti.sh \
  in1="${FWD_READ}" \
  in2="${REV_READ}" \
  sweep=25,31,37,45,55,67,81,91 \
  stdev \
  out="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_kmercountmulti.txt"

bbduk.sh \
  -Xmx${MAX_MEMORY}g \
  t=${NCPUS} \
  in1="${FWD_READ}" \
  in2="${REV_READ}" \
  stats="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_stats.txt" \
  bhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_bhist.txt" \
  qhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_qhist.txt" \
  qchist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_qchist.txt" \
  aqhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_aqhist.txt" \
  bqhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_bqhist.txt" \
  lhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_lhist.txt" \
  gchist="stats/${BASE_NAME}_pair${PAIR_INDEX}_raw_gchist.txt" \
  gcbins="auto"



# STEP 2. Perform adapter trimming and weak quality trimming.

## We're just doing one pair in the array as an example.
## This could be a for-loop or parallelised.
PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"
REV_READ="${R2_FQ[${PAIR_INDEX}]}"

cutadapt \
  --quality-cutoff "${TRIM_PHRED},${TRIM_PHRED}" \
  -a "file:${DATA_DIR}/truseq_fwd.fasta" \
  -A "file:${DATA_DIR}/truseq_rev.fasta" \
  --minimum-length "${MIN_READ_LENGTH}" \
  -n 3 \
  --cores ${NCPUS} \
  -o "tmp1_${FWD_READ}" \
  -p "tmp1_${REV_READ}" \
  "${FWD_READ}" \
  "${REV_READ}" \
  > "stats/${BASE_NAME}_pair${PAIR_INDEX}_cutadapt_pass1.txt"

cutadapt \
  --quality-cutoff "${TRIM_PHRED},${TRIM_PHRED}" \
  -a "file:${DATA_DIR}/truseq_fwd.fasta" \
  -A "file:${DATA_DIR}/truseq_rev.fasta" \
  --minimum-length "${MIN_READ_LENGTH}" \
  -n 3 \
  --cores ${NCPUS} \
  -o "trimmed/${FWD_READ}" \
  -p "trimmed/${REV_READ}" \
  "tmp1_${FWD_READ}" \
  "tmp1_${REV_READ}" \
  > "stats/${BASE_NAME}_pair${PAIR_INDEX}_cutadapt_pass2.txt"

rm -f "tmp1_${FWD_READ}" "tmp1_${REV_READ}"


# Run STEP 1 again using the trimmed reads.
## Make sure to change the output names so you don't overwrite the old ones.


# STEP 3. Filter out potential synthetic contaminants including PHiX.

## We're just doing one pair in the array as an example.
## This could be a for-loop or parallelised.
PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"
REV_READ="${R2_FQ[${PAIR_INDEX}]}"

mkdir filtered_synthetic

bbduk.sh \
  -Xmx${MAX_MEMORY}g \
  t=${NCPUS} \
  in1="${FWD_READ}" \
  in2="${REV_READ}" \
  out1="filtered_synthetic/${FWD_READ}" \
  out2="filtered_synthetic/${REV_READ}" \
  ref="${DATA_DIR}/synthetic_contaminants.fasta" \
  stats="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_stats.txt" \
  bhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_bhist.txt" \
  qhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_qhist.txt" \
  qchist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_qchist.txt" \
  aqhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_aqhist.txt" \
  bqhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_bqhist.txt" \
  lhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_lhist.txt" \
  gchist="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_gchist.txt" \
  gcbins="auto" \
  k=31 \
  hdist=0 \
  mcf=0.7 \
  minavgquality="${FILTER_PHRED}" \
  minlength="${MIN_READ_LENGTH}"

kmercountmulti.sh \
  in1="filtered_synthetic/${FWD_READ}" \
  in2="filtered_synthetic/${REV_READ}" \
  sweep=25,31,37,45,55,67,81,91 \
  stdev \
  out="stats/${BASE_NAME}_pair${PAIR_INDEX}_synthetic_contaminant_filtered_kmercountmulti.txt"


# STEP 4. Check for contaminants using kraken2.

## Kraken has a number of databases it knows how to download from the NCBI.
## I used these ones.
DBS="bacteria archaea protozoa viral UniVec_Core fungi human"

# You can add your own genomes to REFS.
# An example of how to add taxids required by kraken using sed is below.
# Just replace 321614 with the actual ncbi taxid
# sed -r '/^>/{s/(^>[^[:space:]]*)/\1|kraken:taxid|321614/g}' \
#   SN15v9_OM_Chr_and_tigs.fasta \
#   > SN15.fasta
# I used the four reference genomes.
REFS="SN15.fasta SN4.fasta SN2000.fasta SN79.fasta"

# first we build the database

kraken2-build --download-taxonomy --db krakendb

for db in ${DBS}; do
  kraken2-build \
    --threads ${NCPUS} \
    --download-library ${db} \
    --db krakendb
done

for fasta in ${REFS}; do
    kraken2-build \
      --threads ${NCPU} \
      --no-masking \
      --add-to-library ${fasta} \
      --db krakendb
done

kraken2-build \
    --threads ${NCPU} \
    --build \
    --db krakendb \
    --kmer-len 35 \
    --minimizer-len 31 \
    --minimizer-spaces 6

kraken2-build --clean --db krakendb

# Now we can search the database for each set of reads.

## We're just doing one pair in the array as an example.
## This could be a for-loop or parallelised.
PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"
REV_READ="${R2_FQ[${PAIR_INDEX}]}"

kraken2 \
  --threads ${NCPUS} \
  --confidence 0.2 \
  --minimum-base-quality 25 \
  --paired \
  --output "kraken/${BASE_NAME}_pair${PAIR_INDEX}.tsv" \
  --report "kraken/${BASE_NAME}_pair${PAIR_INDEX}_report.txt" \
  --db krakendb \
  "${FWD_READ}" \
  "${REV_READ}"


# If you find a large number of reads assigned something other than your
# organism, you could filter them using bbduk, or exclude the sample.


# STEP 5. Run fastqc on all of the reads.

## We're just doing one pair in the array as an example.
## This could be a for-loop or parallelised.
PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"

for r in "${FWD_READ}" "trimmed/${FWD_READ}" "filtered_synthetic/${FWD_READ}"
do
    fastqc "${r}"
done

# Repeat for the reverse reads

# Move all of the files into the stats folder.
# You could do this with xargs instead of a for loop.
for f in $(find . -name "*_fastqc.html")
    mv "${f}" stats
do

for f in $(find . -name "*_fastqc.zip")
    mv "${f}" stats
do


# STEP 6. Align reads to reference genomes and get stats

# Index the reference
REF="SN15.fasta"

bbmap.sh ref="${REF}"


# Align the reads and get some stats.
PAIR_INDEX=0
FWD_READ="${R1_FQ[${PAIR_INDEX}]}"
REV_READ="${R2_FQ[${PAIR_INDEX}]}"

bbmap.sh \
  -Xmx${MAX_MEMORY}g \
  threads=${NCPUS} \
  in1="${FWD_READ}" \
  in2="${REV_READ}" \
  out="alignments/${BASE_NAME}_pair${PAIR_INDEX}.sam" \
  fast \
  local \
  covstats="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_constats.txt" \
  covhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_covhist.txt" \
  basecov="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_basecov.txt" \
  bincov="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_bincov.txt" \
  bhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_bhist.txt" \
  qhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_qhist.txt" \
  aqhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_aqhist.txt" \
  lhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_lhist.txt" \
  ihist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_ihist.txt" \
  ehist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_ehist.txt" \
  qahist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_qahist.txt" \
  indelhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_indelhist.txt" \
  mhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_mhist.txt" \
  gchist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_gchist.txt" \
  idhist="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_idhist.txt" \
  scafstats="stats/${BASE_NAME}_pair${PAIR_INDEX}_aligned_${REF%%.*}_scafstats.txt" \
  gcbins=auto \
  idbins=auto


# STEP 7. Run multiqc on different combinations of stats.

# Often it's hard to really view everything in one go if you have lots of samples.
# I'd doing 1 for each step (initial stats, cutadapt + stats,
# synthetic contaminant filtered, and aligned to references).
# It's also often useful to plot R1 and R2 separately.

multiqc ./stats --filename "multiqc"

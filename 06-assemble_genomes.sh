#!/usr/bin/env bash

module load java/jdk1.8.0_51
module load bbmap/38.39-bin
module switch PrgEnv-cray/6.0.4 PrgEnv-gnu/6.0.4
module load python/3.6.3
module load java/jdk1.8.0_51
module load quast/5.0.2-bin
module load R/3.5.3
module load spades/3.13.0-bin


bbmerge-auto.sh \
  in1="${IN_DIR}/${ROWS[read1_file]}" \
  in2="${IN_DIR}/${ROWS[read2_file]}" \
  out="${OUT_DIR}/${ROWS[merged_file]}" \
  outu1="${OUT_DIR}/${ROWS[read1_file]}" \
  outu2="${OUT_DIR}/${ROWS[read2_file]}" \
  ihist="${OUT_DIR}/${BASENAME}_insert.txt" \
  vstrict \
  rem \
  k=62 \
  extend2=50 \
  ecct



# function from utils.sh
# Returns associative array "ROWS" containing data
get_table_rows "${READ_FILE}" "${SAMPLE_LINE_NUMBERS}"
# Formats the read files to the parameter format required for spades.
format_read_args () {
    ARR=( $1 )
    READ=$2
    STR=""
    for i in "${!ARR[@]}"; do
        j=$(( i + 1 ))
        STR+=" --pe${j}-${READ} ${IN_DIR}/${ARR[i]}"
    done
    echo ${STR}
    }

# Actually calls the arg formatter.
FWD_ARGS=$(format_read_args "${ROWS[read1_file]}" 1)
REV_ARGS=$(format_read_args "${ROWS[read2_file]}" 2)

echo "CHECK: SPAdes ${SAMPLE} OUT_DIR=${OUT_DIR} FWD_ARGS=${FWD_ARGS} REV_ARGS=${REV_ARGS}"

    --careful -k 31,51,71,81,101,127 \


spades.py \
    -o "${OUT_DIR}/${SAMPLE}" \
    ${FWD_ARGS} \
    ${REV_ARGS} \
    --careful -k 21,31,51,71,81,101 \
    --cov-cutoff auto \
    --memory ${MAX_MEMORY} \
    --threads ${NTHREADS}


# step 3

# function from utils.sh
# Returns associative array "ROWS" containing data
get_table_rows "${READ_FILE}" "\${SAMPLE_LINE_NUMBERS}"
# Formats the read files to the parameter format required for spades.
format_read_args () {
    ARR1=( \$1 )
    ARR2=( \$2 )
    STR=""
    for i in "\${!ARR1[@]}"; do
        j=\$(( i + 1 ))
        STR+=" --pe1 ${FASTQ_DIR}/\${ARR1[i]}"
        STR+=" --pe2 ${FASTQ_DIR}/\${ARR2[i]}"
    done
    echo \${STR}
    }

srun \
  --ntasks=1 \
  --ntasks-per-node=1 \
  --export=all \
  quast.py \
    -R ${GENOME_FILE} \
    -g data/20180327_SN15_Annotation_SetA_TIDY.gff3 \
    -o "${OUT_DIR}/\${SAMPLE}" \
    --min-contig 500 \
    --rna-finding \
    --threads ${NCPU_PER_NODE} \
    --gene-finding \
    --conserved-genes-finding \
    --fungus \
    --labels "\${LABELS}" \
    \$(format_read_args "\${ROWS[read1_file]}" "\${ROWS[read2_file]}") \
    \${SCAFFOLDS[@]} \${CONTIGS[@]}



# function from utils.sh
# Returns associative array "ROWS" containing data
get_table_rows "${READ_FILE}" "${SAMPLE_LINE_NUMBERS}"
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
    -R ${GENOME_FILE} \
    -g data/20180424_SN15_Annotation_SetA_TIDY.gff3 \
    -o "${OUT_DIR}/${SAMPLE}" \
    --min-contig 500 \
    --rna-finding \
    --threads ${NCPU_PER_NODE} \
    --gene-finding \
    --conserved-genes-finding \
    --fungus \
    --labels "${LABELS}" \
    $(format_read_args "${ROWS[read1_file]}" "${ROWS[read2_file]}" "${ROWS[merged_file]}") \
    ${SCAFFOLDS[@]} \${CONTIGS[@]}

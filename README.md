# 20200306-pnodorum_pangenome_pseudocode

Code accompanying an upcoming manuscript on a P. nodorum pangenome.


## Data

In the folder `data/` is a small number of small files used for filtering and adapter trimming.

- `data/truseq_fwd.fasta` and `data/truseq_rev.fasta` both come from the illumina documentation.
- `data/synthetic_contaminants.fasta` contains the PHiX genome (refseq: NC_001422.1) and
   a common truseq PCR dimer product that is included in the bbmap suite.


## `01-quality_control.sh`

This is the initial quality control step for all raw fastq reads.
The outputs is a set of fastq reads that have been trimmed and filtered for synthetic contaminants, and a number of qc statistics.

The pseudocode is based on a [nextflow](https://www.nextflow.io/) pipeline called [qcflow](https://github.com/darcyabjones/qcflow), which is how the commands were actually run.

At a high level, the reads are trimmed using cutadapt, and common contaminants in illumina sequencing (in `data/synthetic_contaminants.fasta`) are filtered out using bbduk.
We also screen for other contamination using kraken2, and align reads to the reference genomes using bbmap to find insert sizes.
Statistics for each step are computed using the bbsuite, fastqc and multiqc.


## `02-align_reads.sh`

This is a script to align reads to a reference genome e.g. SN15.
This uses the trimmed and synthetic contaminant filtered reads from `01-quality_control.sh` as fastq input.
The output is a single BAM file per sample.

We ran a version of this script on the Magnus cluster at [Pawsey](https://pawsey.org.au/), using slurm batch scripts.
The reads are aligned using [stampy](https://www.well.ox.ac.uk/research/research-groups/lunter-group/lunter-group/stampy) separately for each read pair, and then combined and PCR duplicates are marked using [Picard](https://broadinstitute.github.io/picard/).


## `03-predict_variants.sh`

Predict short variants in the genomes using read alignments from `02-align_reads.sh` using something like the GATK bootstrapped procedure.

We don't have a reference set of SNPs to recalibrate the base quality scores with, so we have to perform a first pass without recalibration.
The steps are something like this:

1. Predict haplotypes for each of the samples.
2. Perform join genotyping for all of the samples combined.
3. Compute statistics for each of the variant types separately, and plot them.
4. Filter the vcfs for each of the variant types using hard thresholds determined from the plots in step 3.
5. Recalibrate the base quality scores in the original (unrecalibrated) BAMs using the filtered variants.
   If this is the first pass through the pipeline, go back to step 1, using the recalibrated bam file as input instead of the raw one.
   Otherwise compare the GATK recalibration "covariates" with those from the previous iteration.
   If there is a difference, then go back to step 1 and use the recalibrated BAM file as input instead of the original ones.
   If there is no difference between the covariates, then you can stop the bootstrapping procedure and take the unfiltered vcf from step 3 as the final set (filtered using a different procedure).

This pipeline uses some additional python scripts found in the `bin/` folder.
These require python3 with pandas, matplotlib, and seaborn installed.


## `04-filter_variants.sh`

Perform a final pass of filtering variants with much more in-depth threshold selection and visualisation.

The final genotype thresholds we used were.
MIN_DP=8
MIN_GQ=30

Any variant loci with no non-null genotypes were masked out.

This pipeline uses some additional python and R scripts found in the `bin/` folder.
A conda environment was used to install the R and python packages required, a template for which is available in `./04-filter_variants_environment.yaml`.

To setup the environment install [miniconda](https://docs.conda.io/en/latest/miniconda.html) (or the full anaconda install) and run:

```bash
# This will create an environment called "r_snps"
conda env create -f "04-filter_variants_environment.yaml"
conda activate r_snps

# Run the scripts.

conda deactivate
```


## `05-run_snpeff.sh`

Annotate the variants by their predicted effect on coding genes using SNPEff.

A Dockerfile is provided which sets up SNPEff with P. nodorum SN15 in it's database.

To build the dockerfile, you'd need to install [docker](https://www.docker.com/) (or similar), and have the files `SN15.fasta`, `SN15.faa`, and `SN15.gff3` in the `data/` directory.

Then run:

```bash
sudo docker build -t "darcyabjones/snpeff" --force-rm -f 05-run_snpeff.Dockerfile .
```

Then to run the actual script.

```bash
sudo docker run --rm -v "${PWD}:/data:rw" -w "/data" "darcyabjones/snpeff" 05-run_snpeff.sh
```


## `06-assemble_genomes.sh`

module load java/jdk1.8.0_51
module load bbmap/38.39-bin
module switch PrgEnv-cray/6.0.4 PrgEnv-gnu/6.0.4
module load python/3.6.3
module load java/jdk1.8.0_51
module load quast/5.0.2-bin
module load R/3.5.3
module load spades/3.13.0-bin

Assemble the genomes using spades, and bbmerge.
Also perform qc using quast.


## `07-assemble_mitochondrial_genomes.sh`

TODO: based on mitoflow.

Filter mitochondrial genomes from nuclear genomes (in postasm).

```

darcyabjones/mitoflow
af42920
```

## `08-assembly_qc.sh`

TODO: based on qcflow/postasm


## `09-reference_based_pav.sh`

TODO: based on mumflow.


## `10-transposable_element_prediction.sh`

TODO: based on pante


## `11-gene_prediction.sh`

TODO: based on panann

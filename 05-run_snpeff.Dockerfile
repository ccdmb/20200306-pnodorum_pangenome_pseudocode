FROM quay.io/biocontainers/snpeff:4.3.1t--1

ADD data/SN15.gff3 /usr/local/share/snpeff-4.3.1t-1/data/Pnod.v9/genes.gff
ADD data/SN15.fasta /usr/local/share/snpeff-4.3.1t-1/data/Pnod.v9/sequences.fa
ADD data/SN15.faa /usr/local/share/snpeff-4.3.1t-1/data/Pnod.v9/protein.fa

WORKDIR /usr/local/share/snpeff-4.3.1t-1
RUN set -eu \
    && echo "# Pnod genome, version Pnod.v9" >> snpEff.config \
    && echo "Pnod.v9.genome : Pnod" >> snpEff.config \
    && snpEff build -gff3 -v Pnod.v9


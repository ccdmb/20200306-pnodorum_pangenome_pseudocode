#!/usr/bin/env bash

set -eu

OUTDIR="run_snpeff"
VCF="filtered.vcf.gz"

mkdir "${OUTDIR}"

snpEff -v -nodownload Pnod.v9 "${VCF}" > "${OUTDIR}/$(basename ${VCF})"

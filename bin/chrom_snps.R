#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vcfR"))

option_list <- list(
    make_option(
        c("-v", "--vcf"),
        type="character",
        action="store",
        help="The VCF file to parse (required)."
    ),
    make_option(
        c("-f", "--fasta"),
        type="character",
        action="store",
        help="The fasta file to parse (required)."
    ),
    make_option(
        c("-g", "--gff"),
        type="character",
        action="store",
        help="The gff file to parse (optional)."
    ),
    make_option(
        c("-o", "--outdir"),
        type="character",
        action="store",
        help="Where to store the results (required)."
    ),
    make_option(
        c("-c", "--chromosome"),
        type="character",
        action="store",
        default=NULL,
        help="Select a chromosome to plot. Default will plot each chromosome."
    ),
    make_option(
        "--include-filter",
        dest="include_filter",
        type="character",
        action="store",
        default="PASS",
        help="Value of filter column to retain for plotting."
    ),
    make_option(
        c("-q", "--min-qual"),
        dest="min_qual",
        type="double",
        action="store",
        default=NULL,
        help="Filter qual lower than this [default %default]."
    ),
    make_option(
        c("-d", "--min-dp"),
        dest="min_dp",
        type="double",
        action="store",
        default=NULL,
        help="Filter dp lower than this [default %default]."
    ),
    make_option(
        c("-p", "--max-dp"),
        dest="max_dp",
        type="double",
        action="store",
        default=NULL,
        help="Filter dp greater than this [default %default]."
    ),
    make_option(
        c("-m", "--min-mq"),
        dest="min_mq",
        type="double",
        action="store",
        default=NULL,
        help="Filter mq lower than this [default %default]."
    ),
    make_option(
        c("-u", "--max-mq"),
        dest="max_mq",
        type="double",
        action="store",
        default=NULL,
        help="Filter mq greater than this [default %default]."
    ),
    make_option(
        "--version",
        type="logical",
        action="store_true",
        default=FALSE,
        help="Print extra output.",
    )
)

parser <- OptionParser(
    usage = "%prog --vcf my.vcf --fasta my.fasta --gff my.gff3",
    option_list = option_list
)

args <- parse_args(parser)

get_chr_names <- function(vcf) {
    return(unique(getCHROM(vcf)))
}

get_chr <- function(chromosome, vcf, fasta, gff=NULL) {
    if (isnull(gff)) {
        ann = gff
    } else {
        ann = gff[, gff[1] == chromosome]
    }

    chrom <- create.chromR(
        vcf = vcf[getCHROM(vcf) == chromosome],
        seq = fasta[chromosome],
        ann = ann,
        verbose = FALSE,
    )

    chrom <- proc.chromR(chrom, verbose = FALSE)
    return(chrom)
}

log_stderr <- function(...) {
    cat(sprintf(...), sep='', file=stderr())
}

quit_with_err <- function(...) {
    log_stderr(...)
    quit(save = "no", status = 1, runLast = FALSE)
}

validate_file <- function(path) {
    if (is.null(path)) {
        quit_with_err("Please provide required file")
    }
}

all_filters_null <- function(args) {
    all_null <- is.null(args$min_qual) && is.null(args$min_dp) && is.null(args$max_dp) && is.null(args$min_mq) && is.null(args$max_mq)
    return(all_null)
}

any_filters_null <- function(args) {
    any_null <- is.null(args$min_qual) || is.null(args$min_dp) || is.null(args$max_dp) || is.null(args$min_mq) || is.null(args$max_mq)
    return(any_null)
}



main <- function(args) {
    validate_file(args$vcf)
    validate_file(args$fasta)

    if (is.null(args$outdir)) {
        quit_with_err("Please specify an output directory.")
    }

    # Set defaults for filters.
    if (!all_filters_null(args)) {
        if (is.null(args$min_qual)) {
            args$min_qual = 0
        }

        if (is.null(args$min_dp)) {
            args$min_dp = 0
        }

        if (is.null(args$max_dp)) {
            args$min_dp = Inf
        }

        if (is.null(args$min_mq)) {
            args$min_mq = 0
        }

        if (is.null(args$max_mq)) {
            args$max_mq = Inf
        }
    }

    use_gff = !is.null(args$gff)

    vcf <- read.vcfR(args$vcf)
    vcf <- vcf[getFILTER(vcf) == args$include_filter, ]
    fasta <- ape::read.dna(args$fasta, format = "fasta")

    if (use_gff) {
        gff <- read.table(
            args$gff,
            sep = "\t",
            quote = "",
            stringsAsFactors = FALSE
        )
    }

    available_chrom <- get_chr_names(vcf)
    if (!is.null(args$chromosome)) {
        if ( !(args$chromosome %in% available_chrom) ) {
            quit_with_err("The specified chromosome isn't in the vcf.")
        } else {
            chromosomes <- args$chromosome
        }
    } else {
        chromosomes <- available_chrom
    }


    # Create the output directory
    if (!dir.exists(args$outdir)) {
        dir.create(args$outdir)
    }

    for ( chrom in chromosomes ) {
        if (use_gff) {
            obj <- create.chromR(
                vcf = vcf[getCHROM(vcf) == chrom],
                seq = fasta[chrom],
                ann = gff[gff[1] == chrom, ],
                verbose = FALSE
            )
        } else {
            obj <- create.chromR(
                vcf = vcf[getCHROM(vcf) == chrom],
                seq = fasta[chrom],
                verbose = FALSE
            )
        }

        obj <- proc.chromR(obj, verbose = FALSE)

        pdf(file.path(args$outdir, paste0(chrom, ".pdf")))
        plot(obj)
        dev.off()

        pdf(file.path(args$outdir, paste0(chrom, "_chromoqc.pdf")))
        chromoqc(obj, dp.alpha = 22)
        dev.off()

        if (!all_filters_null(args)) {
            new_obj <- masker(
                obj,
                min_QUAL = args$min_qual,
                min_DP = args$min_dp,
                max_DP = args$max_dp,
                min_MQ = args$min_mq,
                max_MQ = args$max_mq
            )

            new_obj <- proc.chromR(new_obj, verbose = FALSE)

            pdf(file.path(args$outdir, paste0(chrom, "_filtered.pdf")))
            plot(new_obj)
            dev.off()

            pdf(file.path(args$outdir, paste0(chrom, "_filtered_chromoqc.pdf")))
            chromoqc(new_obj, dp.alpha = 22)
            dev.off()
        }
    }
}

main(args)

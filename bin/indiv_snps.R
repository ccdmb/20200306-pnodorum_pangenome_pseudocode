#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("vcfR"))
suppressPackageStartupMessages(library("tibble"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("purrr"))
suppressPackageStartupMessages(library("readr"))

option_list <- list(
  make_option(
    c("-v", "--vcf"),
    type="character",
    action="store",
    help="The VCF file to parse (required)."
    ),
  make_option(
    c("-o", "--outdir"),
    type="character",
    action="store",
    help="Where to store the results (required)."
    ),
  make_option(
    "--no-prefilter",
    dest="no_prefilter",
    type="logical",
    action="store_true",
    default=FALSE,
    help="Don't prefilter the data before computing cutoffs."
    ),
  make_option(
    "--include-filter",
    dest="include_filter",
    type="character",
    action="store",
    default="PASS",
    help="Subset the data to only those with this FILTER value before computing cutoffs [default %default]."
    ),
  make_option(
    "--no-postfilter",
    dest="no_postfilter",
    type="logical",
    action="store_true",
    default=FALSE,
    help="Don't filter the data before plotting the final outputs."
    ),
  make_option(
    c("-d", "--dpq"),
    dest="dp_range",
    type="double",
    action="store",
    default=NULL,
    help="Filter dp by this percentile range size, between 0 and 1. E.G. 0.95 for [0.25, 0.975]. [default %default]."
    ),
  make_option(
    c("-m", "--min-dp"),
    dest="min_dp",
    type="double",
    action="store",
    default=NULL,
    help="Filter dp lower than this [default %default]."
    ),
  make_option(
    c("-g", "--min-gq"),
    dest="min_gq",
    type="double",
    action="store",
    default=NULL,
    help="Filter GQ scores lower than this [default %default]."
    ),
  make_option(
    c("-s", "--gqq"),
    dest="gqq",
    type="double",
    action="store",
    default=NULL,
    help="Filter GQ lower than this percentile, between 0 and 1. [default %default]."
    ),
  make_option(
    "--height",
    dest="height",
    type="double",
    action="store",
    default=10,
    help="Figure height in cm. [default %default]."
    ),
  make_option(
    "--width",
    dest="width",
    type="double",
    action="store",
    default=0.5,
    help="Figure width per sample in cm. [default %default]."
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
  usage = "%prog --vcf my.vcf --outdir results",
  option_list = option_list
)

args <- parse_args(parser)

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
  is.null(args$dpq) && is.null(args$min_dp) && is.null(args$min_gq) && is.null(args$gqq)
}

summarise_stat <- function(tbl, column) {
  tbl %>%
    group_by(sample) %>%
    summarise_at(
      .vars = column,
      .funs = list(
        q001 = partial(quantile, probs=0.001, na.rm = TRUE),
        q005 = partial(quantile, probs=0.005, na.rm = TRUE),
        q01 = partial(quantile, probs=0.01, na.rm = TRUE),
        q025 = partial(quantile, probs=0.025, na.rm = TRUE),
        q05 = partial(quantile, probs=0.05, na.rm = TRUE),
        q10 = partial(quantile, probs=0.10, na.rm = TRUE),
        q25 = partial(quantile, probs=0.25, na.rm = TRUE),
        q50 = partial(quantile, probs=0.50, na.rm = TRUE),
        q75 = partial(quantile, probs=0.75, na.rm = TRUE),
        q90 = partial(quantile, probs=0.90, na.rm = TRUE),
        q95 = partial(quantile, probs=0.95, na.rm = TRUE),
        q975 = partial(quantile, probs=0.975, na.rm = TRUE),
        q99 = partial(quantile, probs=0.99, na.rm = TRUE),
        q995 = partial(quantile, probs=0.995, na.rm = TRUE),
        q999 = partial(quantile, probs=0.999, na.rm = TRUE),
        mean = partial(mean, na.rm = TRUE),
        sd = partial(sd, na.rm = TRUE),
        n = length,
        n_na = function(x) { sum(is.na(x)) }
      )
    )
}

plot_dp <- function(dp) {
  ggplot(dp, aes(x = sample, y = dp)) +
    geom_violin(trim = TRUE, scale = "count", adjust = 1.0) +
    scale_y_continuous(
      trans = scales::log2_trans(),
      breaks = c(1, 10, 100, 1000, 10000),
      minor_breaks = c(1:10, 2:10 * 10, 2:10 * 100, 2:10 * 1000)
      ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}

main <- function(args) {
  validate_file(args$vcf)

  if (is.null(args$outdir)) {
    quit_with_err("Please specify an output directory.")
  }

  if (!is.null(args$dp_range)) {
    if (args$dp_range < 0 || args$dp_range > 1) {
      quit_with_err("DP range is invalid. Must be between 0 and 1.")
    }
  }

  if (!is.null(args$gqq)) {
    if (args$gqq < 0 || args$gqq > 1) {
      quit_with_err("gq percentile is invalid. Must be between 0 and 1.")
    }
  }

  # Create the output directory
  if (!dir.exists(args$outdir)) {
    dir.create(args$outdir)
  }

  vcf <- read.vcfR(args$vcf)
  nsamples <- ncol(vcf@gt) - 1

  if (args$no_prefilter) {
    dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  } else {
    dp <- extract.gt(
      vcf[getFILTER(vcf) == args$include_filter],
      element = "DP",
      as.numeric = TRUE
    )
  }

  dp_tbl <- as_tibble(dp, rownames="pos") %>%
    gather(key = "sample", value = "dp", -pos)

  summarise_stat(dp_tbl, "dp") %>%
    write_tsv(file.path(args$outdir, "dp_summary.tsv"))

  gg <- plot_dp(dp_tbl)
  rm(dp_tbl)

  filename <- file.path(args$outdir, "dp_per_sample.png")
  ggsave(
    filename,
    plot = gg,
    width = args$width * nsamples,
    height = args$height,
    units = "cm"
  )

  if (!is.null(args$dp_range)) {
    q <- (1 - args$dp_range) / 2
    q1 <- q
    q3 <- 1 - q

    dp_quantiles <- apply(
      dp,
      MARGIN=2,
      quantile,
      probs=c(q1, q3),
      na.rm=TRUE
    )
  }
  rm(dp)

  if (args$no_prefilter) {
    gq <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)
  } else {
    gq <- extract.gt(
      vcf[getFILTER(vcf) == args$include_filter],
      element = "GQ",
      as.numeric = TRUE
    )
  }

  filename <- file.path(args$outdir, "gq_per_sample.png")
  png(
    filename,
    width = args$width * nsamples,
    height = args$height,
    units = "cm",
    res = 72
  )
  par( mar = c(8,4,4,2) )
  boxplot(gq, las=2, col=2:5, main="Genotype Quality (GQ)")

  dev.off()

  gq_tbl <- as_tibble(gq, rownames="pos") %>%
    gather(key = "sample", value = "gq", -pos)

  summarise_stat(gq_tbl, "gq") %>%
    write_tsv(file.path(args$outdir, "gq_summary.tsv"))

  rm(gq_tbl)

  if (!is.null(args$gqq)) {

    gq_quantiles <- apply(
      gq,
      MARGIN=2,
      quantile,
      probs=args$gqq,
      na.rm=TRUE
    )
  }

  rm(gq)


  # If none of the filters are set, we can just exit now.
  if (all_filters_null(args)) {
    quit(save = "no", status = 0, runLast = FALSE)
  }

  mask <- is.na(vcf@gt[, -1])
  print(paste("# missing before filter", sum(mask)))

  dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE)

  # Filter by depth quantiles
  if (!is.null(args$dp_range)) {
    # Filter lower bound
    dp_tmp <- sweep(dp, MARGIN=2, FUN = "-", dp_quantiles[1,])
    mask[dp_tmp < 0] <- TRUE
    print(paste("# missing after lower dp quantile", sum(mask)))

    # Filter upper bound
    dp_tmp <- sweep(dp, MARGIN = 2, FUN = "-", dp_quantiles[2,])
    mask[dp_tmp > 0] <- TRUE
    print(paste("# missing after upper dp quantile", sum(mask)))
  }

  # Enforce universal minimum depth.
  mask[dp < args$min_dp] <- TRUE
  rm(dp)
  print(paste("# missing after abs min dp", sum(mask)))

  gq <- extract.gt(vcf, element = "GQ", as.numeric = TRUE)

  # Filter by univesal minimum gq
  if (!is.null(args$min_gq)) {
    mask[gq < args$min_gq] <- TRUE
    print(paste("# missing after abs min gq", sum(mask)))
  }

  if (!is.null(args$gqq)) {
    # Filter lower bound
    gq_tmp <- sweep(gq, MARGIN=2, FUN = "-", gq_quantiles)
    mask[gq_tmp < 0] <- TRUE
    print(paste("# missing after lower gq quantile", sum(mask)))
  }

  rm(gq)

  # Apply the filters, mask is already boolean of correct shape.
  is.na( vcf@gt[, -1][mask] ) <- TRUE
  rm(mask)

  # Show how much we cut out.
  print(vcf)


  # Plot dp information
  if (args$no_postfilter) {
    dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE, mask = TRUE)
  } else {
    dp <- extract.gt(
      vcf[getFILTER(vcf) == args$include_filter],
      element = "DP",
      as.numeric = TRUE,
      mask = TRUE
    )
  }

  gg <- as_tibble(dp, rownames="pos") %>%
    gather(key = "sample", value = "dp", -pos) %>%
    plot_dp()

  filename <- file.path(args$outdir, "dp_per_sample_filtered.png")
  ggsave(
    filename,
    plot = gg,
    width = args$width * nsamples,
    height = args$height,
    units = "cm"
  )
  rm(dp)

  if (args$no_postfilter) {
    gq <- extract.gt(vcf, element = "GQ", as.numeric = TRUE, mask = TRUE)
  } else {
    gq <- extract.gt(
      vcf[getFILTER(vcf) == args$include_filter],
      element = "GQ",
      as.numeric = TRUE,
      mask = TRUE
    )
  }

  filename <- file.path(args$outdir, "gq_per_sample_filtered.png")
  png(
    filename,
    width = args$width * nsamples,
    height = args$height,
    units = "cm",
    res = 72
  )
  par( mar = c(8,4,4,2) )
  boxplot(gq, las=2, col=2:5, main="Genotype Quality (GQ)")

  dev.off()
  rm(gq)


  # Summarise missing data
  missing <- as_tibble(is.na(vcf@gt[, -1]))
  missing$locus <- paste0(getCHROM(vcf), '_', getPOS(vcf))

  gg <- missing %>%
    gather(key = "isolate", value = "missing", -locus) %>%
    group_by(isolate) %>%
    summarise(nmissing = sum(missing)) %>%
    ggplot(aes(x=isolate, y=nmissing)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  filename <- file.path(args$outdir, "missing_genotypes_per_sample.png")
  ggsave(
    filename,
    plot = gg,
    width = args$width * nsamples * 0.75,
    height = args$height,
    units = "cm"
  )

  gg <- missing %>%
    gather(key = "isolate", value = "missing", -locus) %>%
    group_by(locus) %>%
    summarise(nmissing = sum(missing)) %>%
    ggplot(aes(x=locus, y=nmissing)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_blank())

  filename <- file.path(args$outdir, "missing_samples_per_locus.png")
  ggsave(
    filename,
    plot = gg,
    width = args$width * nsamples,
    height = args$height,
    units = "cm"
  )

  rm(missing)
}

main(args)

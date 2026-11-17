#! /usr/bin/env -S Rscript --vanilla

# INPUT ----

args <- commandArgs(trailingOnly = TRUE)
sumstats_file <- args[1]

# SETUP ----

suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("arrow"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("fs"))

# CONVERT & ANNOTATE ----

parquet_file_name <- "formatted_sumstats.parquet"
sumstats <- VariantAnnotation::readVcf(sumstats_file)
# jarl-ignore undesirable_operator: I think this function is OK to use here
sumstats <- MungeSumstats:::vcf2df(
    sumstats,
    add_sample_names = FALSE,
    add_rowranges = FALSE
)
id_split <- tstrsplit(sumstats$ID, r"{:|_|/}")
sumstats$ID <- NULL
sumstats[, `:=`(
    CHR_UCSC = id_split[[1]],
    CHR_ENSEMBL = gsub("^chr", "", id_split[[1]]),
    A1 = id_split[[3]],
    A2 = id_split[[4]],
    BP = id_split[[2]]
)]
sumstats[, `:=`(
    ID_UCSC = paste0(
        CHR_UCSC,
        ":",
        BP,
        "_",
        pmin(id_split[[3]], id_split[[4]]),
        "_",
        pmax(id_split[[3]], id_split[[4]])
    ),
    ID_ENSEMBL = paste0(
        CHR_ENSEMBL,
        ":",
        BP,
        "_",
        pmin(id_split[[3]], id_split[[4]]),
        "_",
        pmax(id_split[[3]], id_split[[4]])
    )
)]
id_cols <- c(
    "CHR_UCSC",
    "CHR_ENSEMBL",
    "BP",
    "END",
    "A1",
    "A2",
    "ID_UCSC",
    "ID_ENSEMBL",
    "SNP"
)
id_cols <- intersect(id_cols, names(sumstats))
remaining_cols <- setdiff(names(sumstats), id_cols)
data.table::setcolorder(sumstats, c(id_cols, remaining_cols))
arrow::write_parquet(
    sumstats,
    parquet_file_name
)

# REFERENCE ----

# https://github.com/Bioconductor/VariantAnnotation/issues/57

#! /usr/bin/env -S Rscript --vanilla

args <- commandArgs(trailingOnly = TRUE)
parquet_file <- args[1]
local_lib <- args[2]

# SETUP ----

suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("fs"))
suppressPackageStartupMessages(library("glue"))
suppressPackageStartupMessages(library("polars", lib.loc = local_lib))

# INPUT ----

message(glue("Processing file: {parquet_file}"))

## HELPER FUNCTIONS ----

pl_calculate_pvalue <- function(beta, se) {
    z <- beta$truediv(se)$abs()
    z$map_batches(
        \(s) {
            as_polars_series(2 * pnorm(as.vector(-s)))
        }
    )
}

pl_calculate_log10_pvalue <- function(beta, se) {
    z <- beta$truediv(se)$abs()
    z$map_batches(
        \(s) {
            as_polars_series(
                -((pnorm(as.vector(-s), log.p = TRUE) + log(2)) / log(10))
            )
        }
    )
}

sumstat <- pl$read_parquet(parquet_file)
file_schema <- sumstat$collect_schema()

## Re-calculate P ----

message("Re-calcuating p-values and adding -log10P column")

sumstat <- sumstat$cast(P = pl$Float64)

# Putting these 2 function calls into one polars call crashes R
sumstat <- sumstat$with_columns(
    P = pl$when(pl$col("P")$eq(0))$then(
        pl_calculate_pvalue(pl$col("BETA"), pl$col("SE"))
    )$otherwise(pl$col("P"))
)

if ("LOG10P" %in% names(file_schema)) {
    sumstat <- sumstat$with_columns(
        LOG10P = pl$when(pl$col("P")$eq(0))$then(
            pl_calculate_log10_pvalue(pl$col("BETA"), pl$col("SE"))
        )$otherwise(pl$col("LOG10P"))
    )
} else {
    sumstat <- sumstat$with_columns(
        LOG10P = pl_calculate_log10_pvalue(pl$col("BETA"), pl$col("SE"))
    )
}

## Remove non-standard chromosomes ----

# check if file was already processed
default_chromosomes_ucsc <- c(
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrM",
    "chrMT"
    # "chrUn" <-not default, will be removed
    # "chrY" <-not default, will be removed
)

default_chromosomes_ensembl <- gsub("^chr", "", default_chromosomes_ucsc)

# save those that were on non-standard chromosomes
sumstats_bad_chr <- sumstat$filter(
    pl$col("CHR_UCSC")$cast(
        pl$Enum(default_chromosomes_ucsc),
        strict = FALSE
    )$is_null()
)

rows_removed <- nrow(sumstats_bad_chr)
if (rows_removed != 0) {
    chr_removed <- sumstats_bad_chr$get_column(
        "CHR_UCSC"
    )$unique()$to_r_vector() |>
        paste(collapse = ", ")
    msg <- glue(
        "Removed n={rows_removed} rows matching chromosome(-s): {chr_removed}"
    )
    message(msg)
} else {
    message("No non-standard chromosomes detected")
}

sumstat <- sumstat$cast(
    CHR_UCSC = pl$Enum(default_chromosomes_ucsc),
    CHR_ENSEMBL = pl$Enum(default_chromosomes_ensembl),
    .strict = FALSE
)$filter(pl$col("CHR_UCSC")$is_null()$not())

## Remove non-integer BP ----

sumstats_bad_bp <- sumstat$filter(
    pl$col("BP")$cast(pl$UInt64, strict = FALSE)$is_null()
)

rows_removed <- nrow(sumstats_bad_bp)
if (rows_removed != 0) {
    msg <- glue(
        "Removed n={rows_removed} non-integer base pair position rows"
    )
    message(msg)
} else {
    message("No non-standard base pair positions detected")
}

sumstat <- sumstat$cast(
    BP = pl$UInt64,
    .strict = FALSE
)$filter(pl$col("BP")$is_null()$not())

## Remove non-standard alleles ----

allele_re <- "^[ACGT]+$"

sumstats_bad_alleles <- sumstat$filter(
    pl$col("A1")$str$contains(allele_re)$not() |
        pl$col("A2")$str$contains(allele_re)$not()
)

rows_removed <- nrow(sumstats_bad_alleles)
if (rows_removed != 0) {
    msg <- glue(
        "Removed n={rows_removed} non-standard alleles"
    )
    message(msg)
} else {
    message("No non-standard alleles detected")
}

sumstat <- sumstat$filter(
    pl$col("A1")$str$contains(allele_re) |
        pl$col("A2")$str$contains(allele_re)
)

unique_a1 <- sumstat$get_column("A1")$unique()$to_r_vector()
unique_a2 <- sumstat$get_column("A2")$unique()$to_r_vector()

sumstat <- sumstat$cast(
    A1 = pl$Enum(unique_a1),
    A2 = pl$Enum(unique_a2)
)

if ("N" %in% names(file_schema)) {
    sumstat <- sumstat$cast(
        N = pl$UInt32
    )
}

# SORT ----

sumstat <- sumstat$sort("CHR_UCSC", "BP")

# SAVE REJECTED PARQUET ----

# combine and save the rejected rows

sumstats_bad_bp <- sumstats_bad_bp$cast(
    CHR_UCSC = pl$String,
    CHR_ENSEMBL = pl$String
)

sumstats_bad_alleles <- sumstats_bad_alleles$cast(
    CHR_UCSC = pl$String,
    CHR_ENSEMBL = pl$String,
    BP = pl$String
)

sumstats_rejected <- pl$concat(
    sumstats_bad_chr,
    sumstats_bad_bp,
    sumstats_bad_alleles
)

rows_removed <- nrow(sumstats_rejected)
if (rows_removed != 0) {
    msg <- glue(
        "Removed total n={rows_removed} rows from the object"
    )
    message(msg)

    n_rejected_sig <- nrow(sumstats_rejected$filter(pl$col("P")$lt(5e-8)))

    if (n_rejected_sig > 0) {
        msg <- glue(
            "Rejected variants include n={n_rejected_sig} with p<5e-8"
        )
        message(msg)
    }

    pheotype_base_dir <- "rejected"
    fs::dir_create(pheotype_base_dir)
    reject_file_name <- "rejected_sumstats.parquet"

    message("Writing out rejected variants")
    sumstats_rejected$write_parquet(
        fs::path(
            pheotype_base_dir,
            reject_file_name
        )
    )
} else {
    message(
        "No rows removed from this object, no rejected variants to save"
    )
}

# SAVE PARQUET ----

message("Writing out new `.parquet` object")
sumstat$write_parquet(parquet_file)

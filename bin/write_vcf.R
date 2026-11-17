#! /usr/bin/env -S Rscript --vanilla

# INPUT ----

args <- commandArgs(trailingOnly = TRUE)
parquet_file <- args[1]
local_lib <- args[2]

# SETUP ----

suppressPackageStartupMessages(
    {
        library(fs)
        library(data.table)
        library("polars", lib.loc = local_lib)
        library(glue)
        library(VariantAnnotation)
    }
)

message(glue("Processing file: {parquet_file}"))
sumstat <- pl$read_parquet(parquet_file)

# SAVE VCF ----

message("Creating `GRanges` object")

# This workflow is adapted from MungeSumstats
gr <- GenomicRanges::makeGRangesFromDataFrame(
    df = as.data.table(sumstat),
    keep.extra.columns = TRUE,
    seqnames.field = "CHR_UCSC",
    start.field = "BP",
    end.field = "BP"
)
suppressMessages(
    suppressWarnings(
        GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    )
)
gr$dummy <- "GWAS"

message("Creating `VRanges` object")

vr <- VariantAnnotation::makeVRangesFromGRanges(
    gr,
    ref.field = "A1",
    alt.field = "A2",
    keep.extra.columns = TRUE,
    sampleNames.field = "dummy"
)

message("Writing new VCF files")

vcf_file <- "formatted_sumstats.vcf"
VariantAnnotation::writeVcf(
    obj = vr,
    filename = vcf_file,
    index = TRUE
)

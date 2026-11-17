#! /usr/bin/env -S Rscript --vanilla

# INPUT ARGS ----

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
    make_option(
        c("--pheno"),
        type = "character",
        help = "Phenotype ID",
        metavar = "PHENOTYPE"
    ),
    make_option(
        c("--build"),
        type = "character",
        help = "Genome build the input file is on",
        metavar = "BUILD"
    ),
    make_option(
        c("--file"),
        type = "character",
        help = "Raw summary statistics file",
        metavar = "FILE"
    ),
    make_option(
        c("--col_headers"),
        type = "character",
        help = "File containing custom column headers",
        metavar = "HEADER_FILE"
    ),
    make_option(
        c("--lib"),
        type = "character",
        help = "Additional R library path",
        metavar = "R_LIB"
    ),
    make_option(
        c("--logp_names"),
        type = "character",
        help = "File containing optional log-p-value column names",
        metavar = "LOGP_Names"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# SETUP ----
suppressPackageStartupMessages(library("MungeSumstats"))
suppressPackageStartupMessages(library("GenomicFiles"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("GenomeInfoDb"))
suppressPackageStartupMessages(library("data.table"))

# INPUT VARIABLES ----
phenotype_name <- opt$pheno
genome_build <- trimws(opt$build)
stopifnot(genome_build %in% c("grch37", "grch38"))

raw_sumstat_file <- opt$file
formatted_sumstats_file <- "formatted_sumstats.vcf.bgz"

custom_sumstatsColHeaders <- opt$col_headers
custom_sumstatsColHeaders <- fread(custom_sumstatsColHeaders)
setDF(custom_sumstatsColHeaders)

snplocs_lib <- opt$lib
.libPaths(c(.libPaths(), snplocs_lib))

log_pvalue_column_names_file <- opt$logp_names
log_pvalue_column_names <- fread(log_pvalue_column_names_file)$col_names

# FORMATTING ----
sumstats <- MungeSumstats::format_sumstats(
    path = raw_sumstat_file,
    save_path = formatted_sumstats_file,
    write_vcf = TRUE,

    rmv_chr = c("Y"),
    chr_style = "UCSC",

    return_data = TRUE,
    return_format = 'vranges',

    tabix_index = FALSE,
    INFO_filter = 0.8,
    sort_coordinates = TRUE,
    impute_beta = TRUE,
    impute_se = TRUE,
    ignore_multi_trait = TRUE,
    convert_small_p = FALSE,

    indels = TRUE,
    drop_indels = FALSE,

    bi_allelic_filter = FALSE,
    flip_frq_as_biallelic = TRUE,

    ref_genome = ifelse(genome_build == "grch37", "GRCh37", "GRCh38"),
    on_ref_genome = TRUE,
    convert_ref_genome = NULL,

    mapping_file = custom_sumstatsColHeaders,
    log_folder = "./",
    log_mungesumstats_msgs = TRUE,
    force_new = TRUE
)

# Output is saved in list when logging is active!
sumstats <- sumstats$sumstats

# Remove all column names that are not expected in a GWAS-VCF
canonical_colnames <- c(
    unique(custom_sumstatsColHeaders$Corrected),
    log_pvalue_column_names
)
non_canonical_colnames <- setdiff(
    names(elementMetadata(sumstats)),
    canonical_colnames
)
if (length(non_canonical_colnames) > 0) {
    for (i in non_canonical_colnames) {
        elementMetadata(sumstats)[[i]] <- NULL
    }
}

# Rename existing Log10P column or create it
format_cols <- names(mcols(sumstats))
lp_old <- format_cols[format_cols %in% log_pvalue_column_names]
lp_new <- "LOG10P"
names(lp_new) <- lp_old

# TODO: In case this is >1, Not sure what to expect or to do
if (length(lp_old) == 1) {
    mcols(sumstats) <- rename(mcols(sumstats), lp_new)

    # Make sure the column is positive (converting log10P to -log10P
    mcols(munged_data)[[lp_new]] <- abs(mcols(munged_data)[[lp_new]])
} else {
    new_lp_col <- list("LOG10P" = -log10(mcols(sumstats)$P))
    mcols(sumstats) <- cbind(
        mcols(sumstats),
        new_lp_col
    )
}

# CLEANUP ----
# This is annoying but there's no way around it.
file.remove(formatted_sumstats_file)

# Save with the correct chromosome coding
GenomeInfoDb::seqlevelsStyle(sumstats) <- "UCSC"

# SAVE ----
VariantAnnotation::writeVcf(
    sumstats,
    paste0("formatted_sumstats.vcf")
)

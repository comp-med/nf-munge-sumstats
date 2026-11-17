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
        c("--cores"),
        type = "integer",
        help = "Number of CPU cores to be used",
        metavar = "N_CPU"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# SETUP ----
suppressPackageStartupMessages(library("MungeSumstats"))
suppressPackageStartupMessages(library("GenomicFiles"))
suppressPackageStartupMessages(library("data.table"))

# INPUT VARIABLES ----
snplocs_lib <- opt$lib
phenotype_name <- opt$pheno
raw_sumstat_file <- opt$file
custom_sumstatsColHeaders <- opt$col_headers
ncpus <- opt$cores

.libPaths(c(.libPaths(), snplocs_lib))

# INPUT PREP ----
custom_sumstatsColHeaders <- fread(custom_sumstatsColHeaders)
setDF(custom_sumstatsColHeaders)

# MAIN ----
file_genome_build <- tryCatch(
    # jarl-ignore undesirable_operator: Why would this be a problem here?
    MungeSumstats:::get_genome_build(
        raw_sumstat_file,
        mapping_file = custom_sumstatsColHeaders,
        header_only = TRUE,
        sampled_snps = 50000,
        dbSNP = 155,
        nThread = as.numeric(ncpus)
    ),
    error = function(cnd) {
        message(
            "Trying to impute genome build by comparing number of ",
            "rsIDs that can be mapped with either GRCh37 or GRCh38 dbSNP."
        )
        # From: `https://support.bioconductor.org/p/9153164/`
        suppressMessages(library(
            SNPlocs.Hsapiens.dbSNP155.GRCh37,
            lib.loc = snplocs_lib
        ))
        suppressMessages(library(
            SNPlocs.Hsapiens.dbSNP155.GRCh38,
            lib.loc = snplocs_lib
        ))
        suppressMessages(library(GenomicRanges))
        dat <- fread(raw_sumstat_file, nrows = 5e4)
        dat <- suppressMessages(
            MungeSumstats::standardise_header(
                dat,
                mapping_file = custom_sumstatsColHeaders,
                return_list = FALSE
            )
        )
        gr <- suppressMessages(MungeSumstats:::to_granges(dat))
        known_snps_b37 <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh37, gr)
        known_snps_b38 <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
        hits_b37 <- findOverlaps(gr, known_snps_b37)
        hits_b38 <- findOverlaps(gr, known_snps_b38)
        if (
            anyDuplicated(queryHits(hits_b37)) ||
                anyDuplicated(queryHits(hits_b38))
        ) {
            warning("some SNPs are mapped to more than 1 known SNP")
        }
        if (length(hits_b37) > length(hits_b38)) {
            message(
                paste0(
                    "Managed to map ",
                    length(hits_b37),
                    "rsIDs to 50k variants with GRCh37 compared to ",
                    length(hits_b38),
                    " for GRCh38. Proceeding with GRCh37."
                )
            )
            return("GRCh37")
        } else if (length(hits_b38) > length(hits_b37)) {
            message(
                paste0(
                    "Managed to map ",
                    length(hits_b38),
                    " rsIDs to 50k variants with GRCh38 compared to ",
                    length(hits_b37),
                    " for GRCh37. Proceeding with GRCh38."
                )
            )
            return("GRCh38")
        } else {
            stop("Could not determine genome build and impute rsIDs!")
        }
    }
)
file_genome_build <- tolower(file_genome_build)
writeLines(
    file_genome_build,
    "genome_build",
    sep = ""
)

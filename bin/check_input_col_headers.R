#! /usr/bin/env -S Rscript --vanilla

# INPUT ARGS ----

suppressPackageStartupMessages(library("optparse"))
option_list <- list(
    make_option(
        c("--corrected_column_names"),
        type = "character"
    ),
    make_option(
        c("--redundant_column_names"),
        type = "character"
    ),
    make_option(
        c("--log_pvalue_column_names"),
        type = "character"
    ),
    make_option(
        c("--column_header_table"),
        type = "character"
    )
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# SETUP ----
library(data.table)
library(MungeSumstats)

corrected_column_names_path <- opt$corrected_column_names
redundant_column_names_path <- opt$redundant_column_names
log_pvalue_column_names_path <- opt$log_pvalue_column_names

# INPUT FILE COLUMN NAMES ----
input_files <- unlist(
    fread(
        opt$column_header_table,
        header = FALSE
    ),
    use.names = FALSE
)
all_colnames <- lapply(input_files, function(file) {
    message(file) # TODO use logger instead
    fread(file)
})
all_colnames <- rbindlist(all_colnames)
all_colnames <- unique(toupper(unlist(all_colnames$colnames)))
length(all_colnames) # TODO: turn into debug msg

# COLUMN NAME MAPPING TABLE ----
unknown_colnames <- setdiff(
    all_colnames,
    MungeSumstats:::sumstatsColHeaders$Uncorrected
)
superfluous_names <- fread(redundant_column_names_path)[[1]]
unknown_colnames <- setdiff(unknown_colnames, superfluous_names)

# CUSTOM MAPPING TABLE ----
matching_table_unmatched <- fread(corrected_column_names_path)

# MAPPING TABLE OVERLAP ----
unknown_colnames <- setdiff(
    unknown_colnames,
    matching_table_unmatched$Uncorrected
)

# CHECK FOR LOG10P VALUES ----

logp_colnames <- fread(log_pvalue_column_names_path, sep = "\t")[[1]]
unknown_colnames <- setdiff(
    unknown_colnames,
    logp_colnames
)

if (length(unknown_colnames) != 0) {
    fwrite(data.table(unknown_colnames), "unknown_colnames.csv")
    message("Unknown column headers detected!")
    message("You might want to add them to either table in `assets`")
    message(paste(unknown_colnames, collapse = ", "))
    stop("Unknown column headers! Manual intervention required")
}

# NEW MAPPING TABLE ----
sumstatsColHeaders <- rbindlist(
    list(
        MungeSumstats:::sumstatsColHeaders,
        matching_table_unmatched
    )
)

fwrite(
    sumstatsColHeaders,
    "sumstatsColHeaders.csv",
    sep = ","
)

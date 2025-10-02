process CHECK_INPUT_COL_HEADERS {
    
    cache 'lenient'
    tag 'singe_execution'
    label 'rProcess'
    
    input:
    path input_file_table
    path r_lib

    output:
    path "sumstatsColHeaders.csv"

    script:
    """
    #! /usr/bin/env Rscript

    # SETUP ----
    library(data.table, lib.loc = "$r_lib")
    library(MungeSumstats, lib.loc = "$r_lib")

    # INPUT FILE COLUMN NAMES ----
    input_files <- unlist(
      fread(
        "$input_file_table",
        header = FALSE
        ),
      use.names = FALSE
    )
    vcf_filter <- !grepl(r"[\\.vcf\$|vcf\\.gz\$]", input_files)
    all_colnames <- lapply(input_files[vcf_filter], function(file) {
      names(fread(file, nrows = 0))
    })
    all_colnames <- toupper(unique(unlist(all_colnames)))
    length(all_colnames)

    # COLUMN NAME MAPPING TABLE ----
    unknown_colnames <- setdiff(
      all_colnames, 
      MungeSumstats:::sumstatsColHeaders\$Uncorrected
    )
    superfluous_names <- c(
      "P.R.", "OR.R.", "Q", "I", "N_BBK", "IS_STRAND_FLIP", "IS_DIFF_AF_GNOMAD", 
      "EFFECTIVE_CASES", "CI_LOWER", "CI_UPPER", "HWEP", "IMPUTATIONINFO", 
      "#KEY", "LOGP", "QUAL", "FILTER", "FORMAT", "LP", "CHISQ_LINREG", 
      "CHISQ_BOLT_LMM_INF", "P_BOLT_LMM_INF", "CHISQ_BOLT_LMM", "END", 
      "HM_CI_LOWER", "HM_CI_UPPER", "HM_CODE", "HM_COORDINATE_CONVERSION",
      "POS_REGION_START_REGION_END"
    )
    unknown_colnames <- setdiff(unknown_colnames, superfluous_names)
     
    # CUSTOM MAPPING TABLE ----
    Uncorrected <- c(
      "N_COHORT", "N_IND", "ALL_META_AF", "INV_VAR_META_BETA", "INV_VAR_META_SEBETA", 
      "INV_VAR_META_P", "INV_VAR_HET_P", "N_DATASET", "CASES", "EFFECTSIZE", 
      "EFFECTSIZE_SD", "NCTRLS", "ISQ_HET", "P_HET", "SI", "A1FREQ", 
      "CHROMSOME", "NCONTROLS", "HM_VARIANT_ID", "HM_RSID", "HM_CHROM", 
      "HM_POS", "HM_OTHER_ALLELE", "HM_EFFECT_ALLELE", "HM_BETA", "HM_ODDS_RATIO", 
      "HM_EFFECT_ALLELE_FREQUENCY", "OR_ALLELE", "NONOR_ALLELE", "OR_SE", "NCONTROLS",
      "POSITION_GRCH37", "EFFECT_ALLELE_FREQ_1000G", "METAL_EFFECT", "METAL_STDERR",
      "METAL_PVALUE", "POOLED_ALT_AF", "FREQ1"
    )
    Corrected <- c(
      "NSTUDY", "N", "FRQ", "BETA", "SE", "P", "HETPVAL", "NSTUDY", 
      "N_CAS", "BETA", "SE", "N_CON", "HETISQT", "HETPVAL", "INFO", 
      "FRQ", "CHR", "N_CON", "SNP", "SNP", "CHR", "BP", "A1", "A2", 
      "BETA", "OR", "FRQ", "A2", "A1", "SE", "N_CON", "BP", "FRQ",
      "BETA", "SE", "P", "FRQ", "FRQ"
    )
    matching_table_unmatched <- data.table(Uncorrected, Corrected)

    # MAPPING TABLE OVERLAP ----
    unknown_colnames <- setdiff(
        unknown_colnames,
        matching_table_unmatched\$Uncorrected
    )
    if (length(unknown_colnames) != 0) {
      message("Unknown column headers detected!")
      message("You might want to add them in the process `CHECK_INPUT_COL_HEADERS`")
      message(paste(unknown_colnames, collapse = ", "))
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
    """

    stub:
    """
    touch sumstatsColHeaders.csv
    """
}

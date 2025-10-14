process GET_INPUT_COL_HEADERS {
    cache 'lenient'
    tag "$file_id"
    label 'rProcess'
    
    input:
    tuple val(file_id), path(input_file_path), path(r_lib)

    output:
    path "column_headers.csv"

    script:
    """
    #! /usr/bin/env Rscript
    r_lib  <- "$r_lib"
    library(data.table, lib.loc = r_lib)

    input_file  <- "$input_file_path"
    file_colnames <- toupper(names(fread(input_file, nrows = 0)))
    fwrite(
        list(
        colnames = file_colnames), 
        "column_headers.csv"
    )
    """

    stub:
    """
    touch column_headers.csv
    """
}

process CHECK_INPUT_COL_HEADERS {
    
    cache 'lenient'
    tag 'singe_execution'
    label 'rProcess'
    
    input:
    path column_header_table
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
        "$column_header_table",
        header = FALSE
        ),
      use.names = FALSE
    )
    all_colnames <- lapply(input_files, function(file) {
      message(file) # TODO use logger instead
      fread(file)
    })
    all_colnames <- rbindlist(all_colnames)
    all_colnames <- unique(toupper(unlist(all_colnames\$colnames)))
    length(all_colnames) # TODO: turn into debug msg

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
      "POS_REGION_START_REGION_END","INV_VAR-HET_P",
      "R2", "Q_PVAL", "GENO", "GENO_A", "GENO_U", "P_HWD", "P_HWD_A", "P_HWD_U",
      "F_MISS",
      "F_MISS_A",
      "F_MISS_U",
      "P_MISS",
      "OVERALL",
      "CHR_X_FREQ_ALLELEB_ALL",
      "CHR_X_FREQ_ALLELEB_MALES",
      "CHR_X_FREQ_ALLELEB_FEMALES",
      "CHR_X_ALL_TOTAL",
      "CHR_X_FREQUENTIST_ADD_PVALUE_MALES",
      "CHR_X_FREQUENTIST_ADD_BETA_1.GENOTYPE.SEX.1_MALES",
      "CHR_X_FREQUENTIST_ADD_SE_1.GENOTYPE.SEX.1_MALES",
      "CHR_X_ALL_TOTAL_MALES",
      "CHR_X_FREQUENTIST_ADD_PVALUE_FEMALES",
      "CHR_X_FREQUENTIST_ADD_BETA_1.GENOTYPE.SEX.2_FEMALES",
      "CHR_X_FREQUENTIST_ADD_SE_1.GENOTYPE.SEX.2_FEMALES",
      "CHR_X_ALL_TOTAL_FEMALES",
      "OR_95L",
      "OR_95U",
      "MULTIANCESTRY_BETA_RAND",
      "MULTIANCESTRY_SE_RAND",
      "MULTIANCESTRY_PVAL_RAND",
      "MULTIANCESTRY_HETQTEST",
      "MULTIANCESTRY_DF_HETQTEST",
      "MULTIANCESTRY_PVAL_HETQTEST",
      "EUROPEAN_ANCESTRY_BETA_RAND",
      "EUROPEAN_ANCESTRY_SE_RAND",
      "EUROPEAN_ANCESTRY_PVAL_RAND",
      "EUROPEAN_ANCESTRY_HETQTEST",
      "EUROPEAN_ANCESTRY_DF_HETQTEST",
      "EUROPEAN_ANCESTRY_PVAL_HETQTEST",
      "SE_GC",
      "OR_L95",
      "OR_U95",
      "P_GC",
      "RSQR",
      "WALDCHISQ",
      "R2_ICOGS",
      "HET_I2",
      "HET_P_VALUE",
      "NAME",
      "UNIQID",
      "ZVAL",
      "NCOHORT",
      "CHISQ_ASSOCIATION",
      "NDF_ASSOCIATION",
      "CHISQ_ANCESTRY_HET",
      "NDF_ANCESTRY_HET",
      "P-VALUE_ANCESTRY_HET",
      "CHISQ_RESIDUAL_HET",
      "NDF_RESIDUAL_HET",
      "P-VALUE_RESIDUAL_HET",
      "LNBF",
      "COMMENTS",
      "META_ANALYSIS",
      "RSQ",
      "EXOME",
      "INFO_UKBB",
      "MODEL",
      "RANGE",
      "NGT",
      "WEIGHT",
      "TEST",
      "NMISS",
      "INFO-SCORE",
      "T",
      "SE_T",
      "ALL_TOTAL",
      "CONTROL_AF",
      "CTRL_FREQ1",
      "P.NI",
      "P.I",
      "LOGOR.SE",
      "P_NOSPA",
      "CONVERGE",
      "IMPUTATION_QUALITY",
      "CONTROLS_MAF",
      "ALL_OR_LOWER",
      "ALL_OR_UPPER",
      "FREQUENTIST_ADD_PVALUE",
      "BAYES_FACTOR",
      "ANNO",
      "ANNOFULL",
      "ERROR",
      "NP",
      "CHISQ_CMH",
      "OR_CMH",
      "L95",
      "U95",
      "CHISQ",
      "F_A",
      "DIRECTION_EFFECTS_COHORTS_ALL",
      "HET_ISQ_ALL",
      "HET_PVALUE_ALL",
      "HETP",
      "CMH P",
      "SNP_ONCO",
      "EFFECT_ALLELE_FREQUENCY_CASES",
      "EFFECT_ALLELE_FREQUENCY_CONTROLS"
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
      "METAL_PVALUE", "POOLED_ALT_AF", "FREQ1",
      "RS_ID_ALL",
      "CASE_AF",
      "NUM_CASES",
      "NUM_CONTROLS",
      "INFO_ALL",
      "ALL_MAF",
      "SNP_POS",
      "#SNP",
      "POSITION_B37",
      "CODED",
      "NON_CODED",
      "CODED_FREQ",
      "RS NUMBER OR CHR:POSITION:[SNP/INDEL]",
      "CASE_FREQ1",
      "EFFECT1",
      "RS_DBSNP147",
      "NEG_LOG_10_P_VALUE",
      "REF_ALLELE",
      "NSAMPLE",
      "BETA_0",
      "SE_0",
      "P-VALUE_ASSOCIATION",
      "NTOTAL",
      "SNP.1",
      "FRQ_A_59851",
      "EFFECT_ALLEL",
      "BASE_PAIR",
      "CASES_MAF",
      "ALL_OR",
      "CUM_EFF_SAMPLE_SIZE",
      "NUMBER_CONTROLS",
      "NUMBER_CASES",
      "VARIANT_ID_HG19",
      "BASE_PAIR_LOCATION_GRCH38",
      "EAF_A1",
      "N_EFFECTIVE_SAMPLESIZE_ALL",
      "EAF_REF",
      "CHR:POSITION", 
      "ORX"
    )
    Corrected <- c(
      "NSTUDY", "N", "FRQ", "BETA", "SE", "P", "HETPVAL", "NSTUDY", 
      "N_CAS", "BETA", "SE", "N_CON", "HETISQT", "HETPVAL", "INFO", 
      "FRQ", "CHR", "N_CON", "SNP", "SNP", "CHR", "BP", "A1", "A2", 
      "BETA", "OR", "FRQ", "A2", "A1", "SE", "N_CON", "BP", "FRQ",
      "BETA", "SE", "P", "FRQ", "FRQ",
      "SNP",
      "FRQ",
      "N_CAS",
      "N_CON",
      "INFO",
      "FRQ",
      "BP",
      "SNP",
      "BP",
      "A2",
      "A1",
      "FRQ",
      "SNP",
      "FRQ",
      "BETA",
      "SNP",
      "P",
      "A2",
      "N",
      "BETA",
      "SE",
      "P",
      "N",
      "SNP",
      "FRQ",
      "A2",
      "BP",
      "FRQ",
      "OR",
      "N",
      "N_CAS",
      "N_CON",
      "SNP",
      "BP",
      "FRQ","N", "FRQ",
      "SNP",
      "OR"
    )
    matching_table_unmatched <- data.table(Uncorrected, Corrected)

    # MAPPING TABLE OVERLAP ----
    unknown_colnames <- setdiff(
        unknown_colnames,
        matching_table_unmatched\$Uncorrected
    )
    if (length(unknown_colnames) != 0) {
      fwrite(data.table(unknown_colnames), "unknown_colnames.csv")
      message("Unknown column headers detected!")
      message("You might want to add them in the process `CHECK_INPUT_COL_HEADERS`")
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
    """

    stub:
    """
    touch sumstatsColHeaders.csv
    """
}

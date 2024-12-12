/*
    * I want to create pre-processed data in both GRCh37 and GRCh38
*/
process GET_GENOME_BUILD {
    
    cache true
    tag "$phenotype_name"
    label 'rProcess'

    input:
    tuple val(phenotype_name), path(raw_sumstat_file), path(r_lib)

    output:
    tuple val(phenotype_name), path("genome_build"), path(raw_sumstat_file)

    script:
    """
    #! /usr/bin/env Rscript

    # SETUP ----
    r_lib = "$r_lib"
    library("MungeSumstats", lib.loc = r_lib)
    library("GenomicFiles", lib.loc = r_lib)

    # INPUT VARIABLES ----
    phenotype_name <- "$phenotype_name"
    raw_sumstat_file <- "$raw_sumstat_file"

    # INPUT PREP ----
    raw_sumstat_file <- list(raw_sumstat_file)
    names(raw_sumstat_file) <- phenotype_name 

    file_genome_build <- get_genome_builds(
        raw_sumstat_file,
        header_only = TRUE,
        sampled_snps = 50000,
        dbSNP = 155,
        nThread = 8 # TODO: Make this adaptable!
    )
    file_genome_build <- tolower(unlist(file_genome_build))
    write.table(
        file_genome_build,
        "genome_build", 
        quote = FALSE, 
        col.names = FALSE,
        row.names = FALSE
    )
    """

    stub:
    """
    """

}

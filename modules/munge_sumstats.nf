// To know whether I want to liftover data or not I need to infer genome build
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
    r_lib <- "$r_lib"
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
    writeLines(
        file_genome_build,
        "genome_build", 
        sep = ""
    )
    """

    stub:
    """
    echo "grch37" > genome_build
    """

}

// Main function that 
process FORMAT_SUMSTATS {
    
    cache true
    tag "$phenotype_name, $genome_build"
    label 'rProcess'

    input:
    tuple
        val(phenotype_name),
        val(genome_build),
        path(raw_sumstat_file), 
        path(custom_col_headers),
        path(r_lib)

    output:
    tuple
        val(phenotype_name),
        val(genome_build),
        path("formatted_sumstats_${genome_build}.vcf.bgz")

    script:
    """
    #! /usr/bin/env Rscript
    
    # SETUP ----
    r_lib <- "$r_lib"
    library("MungeSumstats", lib.loc = r_lib)
    library("GenomicFiles", lib.loc = r_lib)
    library("data.table", lib.loc = r_lib)

    # INPUT VARIABLES ----
    phenotype_name <- "$phenotype_name"
    raw_sumstat_file <- "$raw_sumstat_file"
    genome_build <- "$genome_build"
    custom_sumstatsColHeaders <- "$custom_col_headers"
    """

    stub:
    """
    touch formatted_sumstats_${genome_build}.tsv.gz
    """

}

// process MUNGE_SUMSTATS_LIFTOVER {
// 
// }
// TODO: Write dedicated liftover process using bcftools?


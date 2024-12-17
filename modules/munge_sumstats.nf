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

// Format summary statistics without lifting them over
process FORMAT_SUMSTATS {
    
    cache true
    tag "$phenotype_name, $genome_build"
    label 'rProcess'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        path(raw_sumstat_file), 
        path(custom_col_headers),
        path(r_lib)

    output:
    tuple val(phenotype_name),
        val(genome_build),
        path("formatted_sumstats_${genome_build}.vcf")

    script:
    """
    #! /usr/bin/env Rscript
    
    # SETUP ----
    r_lib <- "$r_lib"
    suppressPackageStartupMessages(library("MungeSumstats", lib.loc = r_lib))
    suppressPackageStartupMessages(library("GenomicFiles", lib.loc = r_lib))
    suppressPackageStartupMessages(library("VariantAnnotation", lib.loc = r_lib))
    suppressPackageStartupMessages(library("GenomeInfoDb", lib.loc = r_lib))
    suppressPackageStartupMessages(library("data.table", lib.loc = r_lib))

    # INPUT VARIABLES ----
    phenotype_name <- "$phenotype_name"
    genome_build <- trimws("$genome_build")
    stopifnot(genome_build %in% c("grch37", "grch38"))

    raw_sumstat_file <- "$raw_sumstat_file"
    formatted_sumstats_file <- paste0(
        "formatted_sumstats_", 
        genome_build,
        ".vcf.bgz"
    )
    custom_sumstatsColHeaders <- "$custom_col_headers"
    custom_sumstatsColHeaders <- fread(custom_sumstatsColHeaders)
    setDF(custom_sumstatsColHeaders)
    dir.create("./logs", showWarnings = FALSE)

    # FORMATTING ----
    sumstats <- MungeSumstats::format_sumstats(
        path = raw_sumstat_file,
        save_path = formatted_sumstats_file,
        write_vcf = TRUE,
        
        rmv_chr = c("Y", "MT"), 
        chr_style = "UCSC",
        
        return_data = TRUE, 
        return_format = 'vranges',
        
        tabix_index = FALSE,
        INFO_filter = 0.8,
        sort_coordinates = TRUE,
        impute_beta = TRUE,
        impute_se = TRUE,
        ignore_multi_trait = TRUE,
        
        indels = TRUE,
        drop_indels = FALSE,
        
        bi_allelic_filter = FALSE,
        flip_frq_as_biallelic = TRUE,
        
        ref_genome = ifelse(genome_build == "grch37", "GRCh37", "GRCH38"), 
        on_ref_genome = TRUE, 
        convert_ref_genome = NULL,
        
        mapping_file = custom_sumstatsColHeaders, 
        nThread = 20,
        log_folder = "./logs/",
        log_mungesumstats_msgs = TRUE,
        imputation_ind = TRUE,
        log_folder_ind = FALSE,
        force_new = TRUE
    )
    
    # Output is saved in list when logging is active!
    sumstats <- sumstats\$sumstats 
    
    # CLEANUP ----
    # This is annoying but there's no way around it.
    file.remove(formatted_sumstats_file)

    # Save with the correct chromosome coding
    GenomeInfoDb::seqlevelsStyle(sumstats) <- "UCSC"

    # SAVE ----
    VariantAnnotation::writeVcf(
        sumstats, 
        paste0("formatted_sumstats_", genome_build, ".vcf")
    )
    """

    stub:
    """
    touch "formatted_sumstats_${genome_build}.vcf"
    """

}

// Indexing and zipping
process SORT_GZIP_INDEX {

    cache true
    tag "$phenotype_name"

    input:
    tuple val(phenotype_name),
        val(genome_build),
        path(formatted_sumstat_file), 
        path(bcftools_liftover_bin),
        path(bgzip_bin)
        
    output:
    tuple val(phenotype_name),
        val(genome_build),
        path("formatted_sumstats_${genome_build}.vcf.gz"),
        path("formatted_sumstats_${genome_build}.vcf.gz.tbi")

    script:
    """
    BCFTOOLS="$bcftools_liftover_bin"
    BGZIP="$bgzip_bin"
    INPUT_VCF="$formatted_sumstat_file"

    # Sort the file (to be sure)
    ./\$BCFTOOLS sort \$INPUT_VCF -o \$INPUT_VCF

    # BGZip the file
    ./\$BGZIP \$INPUT_VCF

    # Index the file
    ./\$BCFTOOLS index --tbi \${INPUT_VCF}.gz
    """

    stub:
    """
    touch "formatted_sumstats_${genome_build}.vcf.gz"
    touch "formatted_sumstats_${genome_build}.vcf.gz.tbi"
    """
  
}

// process GET_LIFTOVER_FILES {
//     
//     cache true
//     tag "single_execution"
// 
//     output:
//     tuple
//         path('hg19.fa'),
//         path('hg38.fa'),
//         path('hg19ToHg38.over.chain.gz'),
//         path('hg38ToHg19.over.chain.gz'),
// 
//     script:
//     """
//     
//     # TODO: make this configuration globally
//     export http_proxy="http://proxy.charite.de:8080"
//     export https_proxy=$http_proxy
//     export HTTPS_PROXY=$http_proxy
//     export HTTP_PROXY=$http_proxy
// 
// 
//     # Chain Files
//     CHAIN1='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
//     CHAIN2='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
// 
//     # Reference Assemblies
//     REF1='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz'
//     REF2='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
// 
//     """
// 
//     stub:
//     """
//     touch hg19ToHg38.over.chain.gz hg38ToHg19.over.chain.gz hg19.fa hg38.fa.gz
//     """
// 
// }


// process LIFTOVER_SUMSTATS {
// 
// }
// TODO: Write dedicated liftover process using bcftools


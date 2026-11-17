process GET_INPUT_COL_HEADERS {
    cache 'lenient'
    tag "$file_id"
    label 'rProcess'
    conda 'conda_envs/r.yml'
    
    input:
    tuple val(file_id), path(input_file_path)

    output:
    path "column_headers.csv"

    script:
    """
    #! /usr/bin/env Rscript --vanilla
    
    library(data.table)
    input_file  <- "$input_file_path"
    file_colnames <- toupper(names(fread(input_file, nrows = 0)))
    fwrite(
        list(
          colnames = file_colnames
        ), 
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
    conda 'conda_envs/r.yml'
    
    input:
    tuple path(column_header_table), 
      path(corrected_column_names), 
      path(redundant_column_names),
      path(log_pvalue_column_names)

    output:
    path "sumstatsColHeaders.csv"

    script:
    """
    check_input_col_headers.R \
      --corrected_column_names $corrected_column_names \
      --redundant_column_names $redundant_column_names \
      --log_pvalue_column_names $log_pvalue_column_names \
      --column_header_table $column_header_table \
    """

    stub:
    """
    touch sumstatsColHeaders.csv
    """
}

// This downloads R packages not available via conda
process DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES {
    cache 'lenient'
    tag "single_execution"
    label 'rProcess'
    conda 'conda_envs/r.yml'
    
    output:
    path "R/"
    path "R/SNPlocs.Hsapiens.dbSNP155.GRCh37"
    path "R/SNPlocs.Hsapiens.dbSNP155.GRCh38"
    path "R/BSgenome.Hsapiens.1000genomes.hs37d5"
    path "R/BSgenome.Hsapiens.NCBI.GRCh38"
    path "R/polars"


    script:
    """
    #! /usr/bin/env Rscript --vanilla

    Sys.setenv(http_proxy  = "$params.proxy_server")
    Sys.setenv(https_proxy = "$params.proxy_server")
    Sys.setenv(HTTPS_PROXY = "$params.proxy_server")
    Sys.setenv(HTTP_PROXY  = "$params.proxy_server")

    library(fs)
    local_lib <- "R"
    fs::dir_create(local_lib)
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh37", lib = local_lib, method = "wget")
    BiocManager::install("SNPlocs.Hsapiens.dbSNP155.GRCh38", lib = local_lib, method = "wget")
    BiocManager::install("BSgenome.Hsapiens.1000genomes.hs37d5", lib = local_lib, method = "wget")
    BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38", lib = local_lib, method = "wget")
    Sys.setenv(NOT_CRAN = "true")
    install.packages("polars", repos = "https://community.r-multiverse.org", lib = local_lib)
    """

    stub:
    """
    mkdir R
    mkdir R/SNPlocs.Hsapiens.dbSNP155.GRCh37
    mkdir R/SNPlocs.Hsapiens.dbSNP155.GRCh38
    mkdir R/BSgenome.Hsapiens.1000genomes.hs37d5
    mkdir R/BSgenome.Hsapiens.NCBI.GRCh38
    mkdir R/polars
    """
}

process DOWNLOAD_OTHER_DATA {

    tag "other: ${phenotype_name}"

    input:
    tuple val(phenotype_name), val(download_link)

    output:
    tuple val(phenotype_name), path('raw_sumstat_file')

    script:
    """
    wget $download_link --output-document raw_sumstat_file
    """

    stub:
    """
    touch raw_sumstat_file
    """

}

process GWAS_CATALOG_SETUP {
    conda 'lftp'
    tag 'gwascatftp'
    label 'rProcess'

    input:
    val r_lib

    output:
    tuple path("harmonized_list"), path("directory_list")

    script:
    """
    #! /usr/bin/env Rscript

    suppressPackageStartupMessages(library(gwascatftp, lib.loc = "$r_lib"))

    gwascat_settings <- gwascatftp::create_lftp_settings(
        lftp_bin = "lftp",
        use_proxy = TRUE, 
        ftp_proxy = "http://proxy.charite.de:8080"
    )
    harmonized_list <- get_harmonised_list(gwascat_settings)
    directory_list <- get_directory_list(gwascat_settings)
    writeLines(harmonized_list, con = "harmonized_list")
    writeLines(directory_list, con = "directory_list")
    """

    stub:
    """
    #! /usr/bin/env Rscript

    suppressPackageStartupMessages(library(gwascatftp, lib.loc = "$r_lib"))

    gwascat_settings <- gwascatftp::create_lftp_settings(
        lftp_bin = "lftp",
        use_proxy = TRUE, 
        ftp_proxy = "http://proxy.charite.de:8080"
    )
    harmonized_list <- get_harmonised_list(gwascat_settings)
    directory_list <- get_directory_list(gwascat_settings)
    writeLines(harmonized_list, con = "harmonized_list")
    writeLines(directory_list, con = "directory_list")
    """

}

process DOWNLOAD_GWAS_CATALOG_DATA {

    tag "gwas_catalog: ${phenotype_name}"

    input:
    tuple val(phenotype_name), val(id)
    each tuple path(harmonized_list), path(directory_list)

    output:
    tuple val(phenotype_name), path('raw_sumstat_file')

    script: 
    """
    #! /usr/bin/env R
    """

    stub:
    """
    touch raw_sumstat_file
    """

}

process DOWNLOAD_OPEN_GWAS_DATA {

    tag "open_gwas: ${phenotype_name}"
    label 'rProcess'

    input:
    tuple val(phenotype_name), val(id)

    output:
    tuple val(phenotype_name), path('raw_sumstat_file')

    script: 
    """
    #! /usr/bin/env R
    """

    stub:
    """
    touch raw_sumstat_file
    """

}



// process PROCESS_NAME {
// 
//   // cache 'lenient'
//   // tag "$input_value"
//   // label 'process_label'
//   // publishDir
// 
//   input:
//   tuple val(some_value), path(some_path)
//   path(some_path)
//   each path(some_path)
// 
//   output:
//   path "./optional_dir/*file.ext", optional: false, emit: output_name
// 
//   script:
//   """
//   ./script.R \
//     $input_parameter \
//     ${task.some_parameter} \
//     ${params.some_parameter}
// 
//   """
// 
// }

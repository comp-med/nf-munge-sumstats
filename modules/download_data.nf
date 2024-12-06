process DOWNLOAD_OTHER_DATA {

    tag "other: ${phenotype_name}"

    input:
    tuple val(phenotype_name), val(download_link)

    output:
    tuple val(phenotype_name), path("raw_sumstat_file.*")

    script:
    """
    # https://stackoverflow.com/questions/965053/extract-filename-and-extension-in-bash
    download_link="$download_link"
    filename=\$(basename -- "\$download_link")
    extension="\${filename##*.}"
    local_name="raw_sumstat_file.\$extension"
    wget \$download_link --output-document \$local_name
    """

    stub:
    """
    touch raw_sumstat_file.gz
    """

}

process GWAS_CATALOG_SETUP {

    conda 'lftp' // TODO create global environment
    tag 'gwascatftp'
    label 'rProcess'

    input:
    path r_lib

    output:
    tuple path("harmonized_list"), path("directory_list")

    script:
    """
    #! /usr/bin/env Rscript

    suppressPackageStartupMessages(library(gwascatftp, lib.loc = "$r_lib"))

    gwascat_settings <- gwascatftp::create_lftp_settings(
        lftp_bin = "lftp",
        use_proxy = TRUE, 
        ftp_proxy = "http://proxy.charite.de:8080" # TODO Make this not depend on charite environment!
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
    label 'rProcess'

    input:
    tuple val(phenotype_name), 
          val(id), 
          path(harmonized_list),
          path(directory_list),
          path(r_lib)

    output:
    tuple val(phenotype_name), path("raw_sumstat_file.*")

    script:
    """
    #! /usr/bin/env Rscript

    suppressPackageStartupMessages(library(data.table, lib.loc = "$r_lib"))
    suppressPackageStartupMessages(library(fs, lib.loc = "$r_lib"))
    suppressPackageStartupMessages(library(gwascatftp, lib.loc = "$r_lib"))

    # load the `gwascatftp` files created in a separate job
    harmonized_list <- unlist(fread("$harmonized_list", header = FALSE))
    directory_list <- unlist(fread("$directory_list", header = FALSE))
    gwas_cat_id <- "$id"

    gwascat_settings <- gwascatftp::create_lftp_settings(
        lftp_bin = "lftp",
        use_proxy = TRUE, 
        ftp_proxy = "http://proxy.charite.de:8080" # TODO make Charite independent
    )

    is_harmonized <- gwascatftp::is_available_harmonised(
      study_accession = gwas_cat_id, harmonized_list
    )
    download_dir <- fs::path("./")
    if (is_harmonized) {
      
      file_link <- gwascatftp::get_harmonised_accession_file_links(
        study_accession = gwas_cat_id, 
        harmonised_list = harmonized_list,
        directory_list = directory_list, 
        list_all_files = FALSE, 
        lftp_settings = gwascat_settings
      )
      
    } else {
      
      file_links <- gwascatftp::get_accession_file_links(
      study_accession = gwas_cat_id, 
      directory_list = directory_list, 
      lftp_settings = gwascat_settings
    )
      # Only keep the .tsv or .tsv.gz
      file_link <- lapply(file_links, function(x){
        x <- x[fs::path_ext(x) %in% c("tsv", "gz")]
        stopifnot(
          "No summary statistic file could be identified for download!" = length(x) > 0
          )
        return(x)
      })
    }

    gwascatftp::download_accession_files_from_ftp(
        accession_file_links = file_link,
        download_directory = download_dir,
        create_accession_directory = FALSE,
        overwrite_existing_files = TRUE, 
        lftp_settings = gwascat_settings
      )

    downloaded_file_name <- fs::path_file(unlist(file_link))
    file_type <- fs::path_ext(downloaded_file_name)
    fs::file_move(
        path = downloaded_file_name, 
        path_ext_set("raw_sumstat_file", file_type)
    )
    """

    stub:
    """
    touch raw_sumstat_file.gz
    """

}

process DOWNLOAD_OPEN_GWAS_DATA {

    tag "open_gwas: ${phenotype_name}"

    input:
    tuple val(phenotype_name), val(id)

    output:
    tuple val(phenotype_name), path('raw_sumstat_file.vcf.gz')

    script: 
    """
    # TODO make this not depend on Charite proxies
    export http_proxy="http://proxy.charite.de:8080"
    export https_proxy=\$http_proxy
    export HTTPS_PROXY=\$http_proxy
    export HTTP_PROXY=\$http_proxy
    wget https://gwas.mrcieu.ac.uk/files/${id}/${id}.vcf.gz \
        --output-document raw_sumstat_file.vcf.gz
    """

    stub:
    """
    touch raw_sumstat_file.vcf.gz
    """

}

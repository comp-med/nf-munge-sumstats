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
    conda 'conda_envs/r.yml'
    
    input:
    tuple path(column_header_table), 
      path(corrected_column_names), 
      path(redundant_column_names)

    output:
    path "sumstatsColHeaders.csv"

    script:
    """
    #! /usr/bin/env Rscript --vanilla

    # SETUP ----
    library(data.table)
    library(MungeSumstats)

    corrected_column_names_path <- "$corrected_column_names"
    redundant_column_names_path <- "$redundant_column_names"

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
    superfluous_names <- fread(redundant_column_names_path)[[1]]
    unknown_colnames <- setdiff(unknown_colnames, superfluous_names)
     
    # CUSTOM MAPPING TABLE ----
    matching_table_unmatched <- fread(corrected_column_names_path)

    # MAPPING TABLE OVERLAP ----
    unknown_colnames <- setdiff(
        unknown_colnames,
        matching_table_unmatched\$Uncorrected
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
    """

    stub:
    """
    mkdir R
    mkdir R/SNPlocs.Hsapiens.dbSNP155.GRCh37
    mkdir R/SNPlocs.Hsapiens.dbSNP155.GRCh38
    mkdir R/BSgenome.Hsapiens.1000genomes.hs37d5
    mkdir R/BSgenome.Hsapiens.NCBI.GRCh38
    """
}

// This modifies the bcftools conda environment to contain the liftover plugin
process SETUP_LIFTOVER {
    tag 'single_execution'
    conda 'conda_envs/bcftools_liftover.yml'

    output:
    file('bcftools/plugins/liftover.so')

    script:
    """
    BCFTOOLS_VER=\$(bcftools --version | head -1 | awk '{print \$2}')

    # Download bcftools and the plugin
    wget https://github.com/samtools/bcftools/releases/download/\${BCFTOOLS_VER}/bcftools-\${BCFTOOLS_VER}.tar.bz2
    tar xjf bcftools-\${BCFTOOLS_VER}.tar.bz2
    mv bcftools-\${BCFTOOLS_VER} bcftools
    cd bcftools
    wget -P plugins https://raw.githubusercontent.com/freeseek/score/master/liftover.c

    # Give defaults to avoid unbound variable error
    : "\${CPATH=}"
    : "\${LIBRARY_PATH=}"
    : "\${LD_LIBRARY_PATH=}"

    # Set env for compilation
    export CPATH=\$CONDA_PREFIX/include:\$CPATH
    export LIBRARY_PATH=\$CONDA_PREFIX/lib:\$LIBRARY_PATH
    export LD_LIBRARY_PATH=\$CONDA_PREFIX/lib:\$LD_LIBRARY_PATH

    # Compile plugin
    make plugins/liftover.so

    mkdir -p \$CONDA_PREFIX/libexec/bcftools
    cp plugins/liftover.so \$CONDA_PREFIX/libexec/bcftools/
    export BCFTOOLS_PLUGINS=\$CONDA_PREFIX/libexec/bcftools
    mkdir -p \$CONDA_PREFIX/etc/conda/activate.d
    echo "export BCFTOOLS_PLUGINS=\$CONDA_PREFIX/libexec/bcftools" >> \$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh
    """

    stub:
    """
    mkdir -p bcftools/plugins
    touch bcftools/plugins/liftover.so
    """
}

process CHECK_LIFTOVER_PLUGIN {
    tag 'single_execution'
    conda 'conda_envs/bcftools_liftover.yml'

    input:
    file(plugin_file)

    output:
    stdout emit: plugin_check

    script:
    """
    bcftools plugin -lv 2>&1 | \
        grep -q liftover && echo "liftover plugin OK" || \
        (echo "liftover plugin MISSING" && exit 1)
    """

    stub:
    """
    echo "liftover plugin OK"
    """
}

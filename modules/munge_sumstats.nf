// To know whether I want to liftover data or not I need to infer genome build
process GET_GENOME_BUILD {
    
    cache 'lenient'
    tag "$phenotype_name"
    label 'rProcess'

    input:
    tuple val(phenotype_name), path(raw_sumstat_file), path(custom_col_headers), path(r_lib)

    output:
    tuple val(phenotype_name), path("genome_build"), path(raw_sumstat_file)

    script:
    """
    #! /usr/bin/env Rscript

    # SETUP ----
    r_lib <- "$r_lib"
    suppressPackageStartupMessages(library("MungeSumstats", lib.loc = r_lib))
    suppressPackageStartupMessages(library("GenomicFiles", lib.loc = r_lib))
    suppressPackageStartupMessages(library("data.table", lib.loc = r_lib))

    # INPUT VARIABLES ----
    phenotype_name <- "$phenotype_name"
    raw_sumstat_file <- "$raw_sumstat_file"
    custom_sumstatsColHeaders <- "$custom_col_headers"
    
    # INPUT PREP ----
    custom_sumstatsColHeaders <- fread(custom_sumstatsColHeaders)
    setDF(custom_sumstatsColHeaders)
    
    # MAIN ----
    file_genome_build <- tryCatch(
        MungeSumstats:::get_genome_build(
            raw_sumstat_file,
            mapping_file = custom_sumstatsColHeaders,
            header_only = TRUE,
            sampled_snps = 50000,
            dbSNP = 155,
            nThread = as.numeric("$task.cpus")
        ),
        error = function(cnd) {
            message(
              "Trying to impute genome build by comparing number of ",
              "rsIDs that can be mapped with either GRCh37 or GRCh38 dbSNP."
            )
            # From: `https://support.bioconductor.org/p/9153164/`
            suppressMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh37))
            suppressMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh38))
            suppressMessages(library(GenomicRanges))
            dat <- fread(raw_sumstat_file, nrows = 5e4)
            dat <- suppressMessages(
              MungeSumstats::standardise_header(
                dat,
                mapping_file = custom_sumstatsColHeaders, 
                return_list = FALSE
              ))
            gr <- suppressMessages(MungeSumstats:::to_granges(dat))
            known_snps_b37 <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh37, gr)
            known_snps_b38 <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)
            hits_b37 <- findOverlaps(gr, known_snps_b37)
            hits_b38 <- findOverlaps(gr, known_snps_b38)
            if (
              anyDuplicated(queryHits(hits_b37)) || 
              anyDuplicated(queryHits(hits_b38))
            ) {
              warning("some SNPs are mapped to more than 1 known SNP")
            }
            if (length(hits_b37) > length(hits_b38)) {
              message(
                paste0(
                  "Managed to map ", 
                  length(hits_b37), 
                  "rsIDs to 50k variants with GRCh37 compared to ", 
                  length(hits_b38), 
                  " for GRCh38. Proceeding with GRCh37."
                )
              )
              return("GRCh37")
            } else if (length(hits_b38) > length(hits_b37)) {
              message(
                paste0(
                  "Managed to map ", 
                  length(hits_b38), 
                  " rsIDs to 50k variants with GRCh38 compared to ", 
                  length(hits_b37), 
                  " for GRCh37. Proceeding with GRCh38."
                )
              )
              return("GRCh38") 
            } else {
              stop("Could not determine genome build and impute rsIDs!")
            }
        }
    )
    file_genome_build <- tolower(file_genome_build)
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
    
    cache 'lenient'
    tag "$phenotype_name, $genome_build"
    label 'rProcess'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path(raw_sumstat_file), 
        path(custom_col_headers),
        path(r_lib)

    output:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
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
        
        ref_genome = ifelse(genome_build == "grch37", "GRCh37", "GRCh38"), 
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
    
    # Remove all column names that are not expected in a GWAS-VCF
    canonical_colnames <- unique(custom_sumstatsColHeaders\$Corrected)
    non_canonical_colnames <- setdiff(
      names(elementMetadata(sumstats)),
      canonical_colnames
    )
    if (length(non_canonical_colnames) > 0) {
      for (i in non_canonical_colnames)
          elementMetadata(sumstats)[[i]] <- NULL
    }

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

// Get chain files and reference sequences necessary for liftover
process GET_LIFTOVER_FILES {
    
    cache 'lenient'
    tag 'single_execution'

    input:
    path(bgzip_bin)

    output:
    tuple path('hg19.fa'),
        path('hg38.fa'),
        path('hg19ToHg38.over.chain.gz'),
        path('hg38ToHg19.over.chain.gz')

    script:
    """
    
    # TODO: make this configuration globally
    export http_proxy="http://proxy.charite.de:8080"
    export https_proxy=\$http_proxy
    export HTTPS_PROXY=\$http_proxy
    export HTTP_PROXY=\$http_proxy

    # Chain Files
    CHAIN1='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
    CHAIN2='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
    wget \$CHAIN1
    wget \$CHAIN2

    # Reference Assemblies
    REF1='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz'
    REF2='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
    wget \$REF1
    wget \$REF2

    # Unzip the reference sequences and remove the the 
    ./bgzip -d  \$(echo "\${REF1##*/}")
    ./bgzip -d  \$(echo "\${REF2##*/}")
    """

    stub:
    """
    touch hg19ToHg38.over.chain.gz hg38ToHg19.over.chain.gz hg19.fa hg38.fa
    """

}

// Liftover summary statistics from given assembly to the other one
process LIFTOVER_SUMSTATS {

    cache 'lenient'
    tag "$phenotype_name:$genome_build->$other_genome_build"
    publishDir = [
        [
            path: { "${params.outDir}/formatted/${phenotype_name}/grch37/" },
            mode: 'copy',
            pattern: "formatted_sumstats_grch37.vcf{.gz,.gz.tbi}"
        ],
        [
            path: { "${params.outDir}/formatted/${phenotype_name}/grch38/" },
            mode: 'copy',
            pattern: "formatted_sumstats_grch38.vcf{.gz,.gz.tbi}"
        ]
    ]

    input:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path(formatted_sumstats),
        path(hg19_reference),
        path(hg38_reference),
        path(hg19_to_38_chain_file),
        path(hg38_to_19_chain_file),
        path(bcftools_liftover_bin),
        path(bgzip_bin)

    output:
    tuple val(phenotype_name),
        path("formatted_sumstats_grch37.vcf.gz"),
        path("formatted_sumstats_grch37.vcf.gz.tbi"),
        path("formatted_sumstats_grch38.vcf.gz"),
        path("formatted_sumstats_grch38.vcf.gz.tbi")

    script:
    """
    # INPUT FILES ----
    HG19_REF="$hg19_reference"
    HG38_REF="$hg38_reference"
    HG19_TO_HG38_CHAIN="$hg19_to_38_chain_file"
    HG38_TO_HG19_CHAIN="$hg38_to_19_chain_file"
    FROM_GENOME_BUILD="$genome_build"
    TO_GENOME_BUILD="$other_genome_build"
    INPUT_VCF="$formatted_sumstats"
    OUTPUT_VCF="formatted_sumstats_\${TO_GENOME_BUILD}.vcf.gz"

    # SORT & BGZIP INPUT ---- 

    # Sort the file (to be sure)
    ./bcftools sort \$INPUT_VCF -o \$INPUT_VCF

    # BGZip the file
    ./bgzip \$INPUT_VCF
    INPUT_VCF="\${INPUT_VCF}.gz"

    # LIFTOVER DIRECTION ----

    if [ "\$FROM_GENOME_BUILD" == "grch37" ]; then
        CHAIN=\$HG19_TO_HG38_CHAIN
        SOURCE_REF=\$HG19_REF
        TARGET_REF=\$HG38_REF
    else
        CHAIN=\$HG38_TO_HG19_CHAIN
        SOURCE_REF=\$HG38_REF
        TARGET_REF=\$HG19_REF    
    fi

    # COLLAPSE VCF ----

    # Create a collapsed VCF and check REF/ALT alignment in the process
    ./bcftools norm --no-version -Ou -m+ \$INPUT_VCF \
    --check-ref ws \
    -f \$SOURCE_REF \
    -o input_COLLAPSED.vcf.gz

    # FIX VCF REF ----

    # Check Allele mismatches etc
    ./bcftools +fixref input_COLLAPSED.vcf.gz -- -f \$SOURCE_REF

    # LIFTOVER ----

    # Liftover
    ./bcftools +liftover --no-version input_COLLAPSED.vcf.gz -Ou -o output_COLLAPSED.vcf.gz -- \
    -s \$SOURCE_REF \
    -f \$TARGET_REF \
    -c \$CHAIN  \
    --reject reject.vcf.bgz \
    --reject-type z \
    --write-src

    # CHECKS ----

    # Check the file using bcftools
    ./bcftools +af-dist output_COLLAPSED.vcf.gz
    ./bcftools +fixref output_COLLAPSED.vcf.gz -- -f \$TARGET_REF
    ./bcftools stats output_COLLAPSED.vcf.gz

    # Sort, expand, index and save
    ./bcftools  sort -Oz output_COLLAPSED.vcf.gz | \
        ./bcftools  norm --no-version -Oz -m- -o output.vcf.gz

    # INFO Column of output file contains content that his difficult
    # to parse downstream. I'll remove it
    ./bcftools annotate -x INFO output.vcf.gz -Oz -o \$OUTPUT_VCF

    # This is because bcftools does not modify in-place, apparently.
    # Clean up
    rm output.vcf.gz
    rm output_COLLAPSED.vcf.gz
    rm input_COLLAPSED.vcf.gz

    # create index
    ./bcftools index -f --tbi \$INPUT_VCF
    ./bcftools index -f --tbi \$OUTPUT_VCF
    """
    
    stub:
    """
    touch "formatted_sumstats_grch37.vcf.gz"
    touch "formatted_sumstats_grch37.vcf.gz.tbi"
    touch "formatted_sumstats_grch38.vcf.gz"
    touch "formatted_sumstats_grch38.vcf.gz.tbi"
    """

}

// To make downstream processing easier, save as tabular parquet
process SAVE_PARQUET {

    cache 'lenient'
    tag "$phenotype_name"
    label 'rProcess'
    publishDir = [
        [
            path: { "${params.outDir}/formatted/${phenotype_name}/grch37/" },
            mode: 'copy',
            pattern: "formatted_sumstats_grch37.parquet"
        ],
        [
            path: { "${params.outDir}/formatted/${phenotype_name}/grch38/" },
            mode: 'copy',
            pattern: "formatted_sumstats_grch38.parquet"
        ]
    ]

    input:
    tuple val(phenotype_name),
        path(formatted_sumstats_grch37),
        path(formatted_sumstats_grch37_index),
        path(formatted_sumstats_grch38),
        path(formatted_sumstats_grch38_index),
        path(r_lib)

    output:
    tuple val(phenotype_name),
        path("formatted_sumstats_grch37.parquet"),
        path("formatted_sumstats_grch38.parquet")


    script:
    """
    #! /usr/bin/env Rscript

    # SETUP ----

    r_lib <- "$r_lib"
    suppressPackageStartupMessages(library("data.table", lib.loc = r_lib))
    suppressPackageStartupMessages(library("arrow", lib.loc = r_lib))
    suppressPackageStartupMessages(library("VariantAnnotation", lib.loc = r_lib))
    suppressPackageStartupMessages(library("fs", lib.loc = r_lib))

    # INPUT ----

    sumstats_grch37_file <- "$formatted_sumstats_grch37"
    sumstats_grch38_file <- "$formatted_sumstats_grch38"

    # CONVERT & ANNOTATE ----

    sumstats_files <- c(
      sumstats_grch37_file, 
      sumstats_grch38_file
    )

    for (i in seq_along(sumstats_files)) {
      
      parquet_file_name <- fs::path_ext_set(
        fs::path_ext_remove(
          fs::path_ext_remove(
            fs::path_file(
              sumstats_files[i]
            )
          )
        ), "parquet"
      )
      
    sumstats <- VariantAnnotation::readVcf(sumstats_files[i])
    sumstats <- MungeSumstats:::vcf2df(
      sumstats, 
      add_sample_names = FALSE, 
      add_rowranges = FALSE
    )
    id_split <- tstrsplit(sumstats\$ID, r"{:|_|/}")
    sumstats\$ID <- NULL
    sumstats[, `:=`(
      CHR_UCSC = id_split[[1]],
      CHR_ENSEMBL = gsub("^chr", "", id_split[[1]]),
      A1 = id_split[[3]],
      A2 = id_split[[4]],
      BP = id_split[[2]]
    )]
    sumstats[, `:=`(
      ID_UCSC = paste0(
        CHR_UCSC,
        ":",
        BP,
        "_",
        pmin(id_split[[3]], id_split[[4]]), 
        "_",
        pmax(id_split[[3]], id_split[[4]])
      ),
      ID_ENSEMBL = paste0(
        CHR_ENSEMBL, 
        ":",
        BP,
        "_",
        pmin(id_split[[3]], id_split[[4]]),
        "_",
        pmax(id_split[[3]], id_split[[4]])
      )
    )]
    id_cols <- c(
      "CHR_UCSC", "CHR_ENSEMBL", "BP", "END", "A1",
      "A2", "ID_UCSC", "ID_ENSEMBL", "SNP"
    )
    id_cols <- intersect(id_cols, names(sumstats))
    remaining_cols <- setdiff(names(sumstats), id_cols)
    data.table::setcolorder(sumstats, c(id_cols, remaining_cols))  
    arrow::write_parquet(
        sumstats,
        parquet_file_name, 
      )
    }
    # REFERENCE ----

    # https://github.com/Bioconductor/VariantAnnotation/issues/57
    """

    stub:
    """
    touch formatted_sumstats_grch37.parquet formatted_sumstats_grch38.parquet
    """

}

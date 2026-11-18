include {
    GET_GENOME_BUILD;
    FORMAT_SUMSTATS;
    GET_LIFTOVER_FILES;
    LIFTOVER_SUMSTATS;
    NO_LIFTOVER_SUMSTATS;
    CREATE_PARQUET;
    CHECK_PARQUET;
    SAVE_VCF
} from '../modules/munge_sumstats.nf'

workflow MUNGE_SUMSTATS {

    take: 
    input_files_ch
    custom_col_headers
    snplocs_lib
    log_pvalue_column_names

    main:

    input_files_ch = GET_GENOME_BUILD (
        input_files_ch
            .combine(custom_col_headers)
            .combine(snplocs_lib)
    ).map {
        tup -> [
        tup[0],
        file(tup[1]).text.replaceAll("\\s",""), 
        tup[2]
    ]}

    // Define other genome build to be used later for liftover
    input_files_ch = input_files_ch.map {
        tup -> [
            tup[0],
            tup[1],
            tup[1] == "grch37" ? "grch38" : "grch37",
            tup[2]
        ]
    }

    // Format files in their current genome build
    FORMAT_SUMSTATS (
        input_files_ch
            .combine(custom_col_headers)
            .combine(snplocs_lib)
            .combine(log_pvalue_column_names)
    )
    formatted_files_ch = FORMAT_SUMSTATS.out.files

    // Make sure any genome build is selected
    if (params.publish_grch37 == false & params.publish_grch38 == false) {
        log.error("Need select at least one genome build for output")
        exit 1

    // When data should be published as grch38, do not run
    // liftover grch38 -> grch37
    } else if (params.publish_grch37 == false)  {
        formatted_files_ch = formatted_files_ch.branch {
            tup -> 
                liftover: tup[1] == "grch37" 
                no_liftover: tup[1] == "grch38" 
        }

    // Reverse: liftover grch38 files, keep grch37 files
    } else if (params.publish_grch38 == false)  {
        formatted_files_ch = formatted_files_ch.branch {
            tup -> 
                liftover: tup[1] == "grch38"
                no_liftover: tup[1] == "grch37"
        }

    // when both formats are selected, liftover everything
    } else {
        formatted_files_ch = formatted_files_ch.branch {
            _tup -> 
                liftover: true
                no_liftover: false 
        }
    }

    // TODO: Only do this when there are files to liftover
    // Liftover required chain files and reference sequences. 
    // These are downloaded in this process
    liftover_files_ch = GET_LIFTOVER_FILES ()

    if (params.publish_vcf == false & params.publish_parquet == false) {
        log.error("Please select at least one output format in")
        exit 1
    }
    
    // Liftover formatted summary statistics from one genome build into the other
    // so either grch37 -> grch38 or the other way around
    LIFTOVER_SUMSTATS (
        formatted_files_ch.liftover
            .combine(liftover_files_ch)
    )
    
    liftover_ch = LIFTOVER_SUMSTATS.out.files.flatMap { 
        phenotype_name, 
        _genome_build, 
        _other_genome_build, 
        grch37_vcf, 
        grch38_vcf -> 
        [ 
          [phenotype_name, "grch37", grch37_vcf], 
          [phenotype_name, "grch38", grch38_vcf] 
        ]
    }

    // remove files I do not want to publish
    if (params.publish_grch37 == false) {
        liftover_ch = liftover_ch.filter {
            t -> t[1] != "grch37"
        }
    } else if (params.publish_grch38 == false) {
        liftover_ch = liftover_ch.filter {
            t -> t[1] != "grch38"
        }
    }

    // some bcftools checks
    no_liftover_ch = NO_LIFTOVER_SUMSTATS(
        formatted_files_ch.no_liftover.combine(liftover_files_ch)
    )

    liftover_log_ch = LIFTOVER_SUMSTATS.out.log.concat(no_liftover_ch.log)

    // combine the channels again
    liftover_ch = liftover_ch.concat(no_liftover_ch.files)
    
    // Having the files available in a tabular format makes it more convenient
    // Therefore, saving everything in parquet makes sense
    CREATE_PARQUET (liftover_ch)
    CHECK_PARQUET ( CREATE_PARQUET.out.combine(snplocs_lib) )
    SAVE_VCF( CHECK_PARQUET.out.parquet.combine(snplocs_lib) )

    emit:
    parquet = CHECK_PARQUET.out.parquet
    check_rejected = CHECK_PARQUET.out.rejected
    vcf = SAVE_VCF.out
    format_log = FORMAT_SUMSTATS.out.log
    liftover_log = liftover_log_ch
    liftover_rejected = LIFTOVER_SUMSTATS.out.reject

}

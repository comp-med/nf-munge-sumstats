include {
    GET_GENOME_BUILD;
    FORMAT_SUMSTATS;
    GET_LIFTOVER_FILES;
    LIFTOVER_SUMSTATS;
    NO_LIFTOVER_SUMSTATS;
    SAVE_PARQUET;
} from '../modules/munge_sumstats.nf'

workflow MUNGE_SUMSTATS {

    take: 
    input_files_ch
    custom_col_headers
    snplocs_lib

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
    formatted_files_ch = FORMAT_SUMSTATS (
        input_files_ch
            .combine(custom_col_headers)
            .combine(snplocs_lib)
    )

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
    liftover_ch = LIFTOVER_SUMSTATS (
        formatted_files_ch.liftover
            .combine(liftover_files_ch)
    ).flatMap { phenotype_name, 
                _genome_build, 
                _other_genome_build, 
                grch37_vcf, 
                grch37_index, 
                grch38_vcf, 
                grch38_index, 
                _log -> 
        [ 
          [phenotype_name, "grch37", grch37_vcf, grch37_index], 
          [phenotype_name, "grch38", grch38_vcf, grch38_index] 
        ]
        }

    // some bcftools checks
    no_liftover_ch = NO_LIFTOVER_SUMSTATS(
        formatted_files_ch.no_liftover
    )

    // Having the files available in a tabular format makes it more convenient
    // Therefore, saving everything in parquet makes sense
    if (params.publish_parquet == true) {

        // combine the channels again
        liftover_ch = liftover_ch.concat(no_liftover_ch)

        SAVE_PARQUET (liftover_ch)
    }
}

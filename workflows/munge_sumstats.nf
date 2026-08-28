include {
    GET_GENOME_BUILD;
    FORMAT_SUMSTATS;
    GET_LIFTOVER_FILES;
    LIFTOVER_SUMSTATS;
    SAVE_PARQUET;
} from '../modules/munge_sumstats.nf'

workflow MUNGE_SUMSTATS {

    take: 
    input_files_ch
    custom_col_headers
    snplocs_lib
    bcftools_liftover_bin

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
    )

    // Liftover required chain files and reference sequences. 
    // These are downloaded in this process
    liftover_files_ch = GET_LIFTOVER_FILES (
    )
    
    // Liftover formatted summary statistics from one genome build into the other
    // so either grch37 -> grch38 or the other way around
    liftover_ch = LIFTOVER_SUMSTATS (
        formatted_files_ch
            .combine(liftover_files_ch)
            .combine(bcftools_liftover_bin)
    )
    
    // Having the files available in a tabular format makes it more convenient
    // Therefore, saving everything in parquet makes sense
    SAVE_PARQUET (
        liftover_ch
    )

    // emit:
    // formatted_data = liftover_ch.out
}

include {
    GET_GENOME_BUILD;
    FORMAT_SUMSTATS
} from '../modules/munge_sumstats.nf'

workflow MUNGE_SUMSTATS {

    take: 
    input_files_ch
    custom_col_headers
    r_lib

    main:

    input_files_ch = GET_GENOME_BUILD (
        input_files_ch.combine(r_lib) 
    ).map {
        tup -> [tup[0], file(tup[1]).text, tup[2]]
    }

    // TODO
    // Run once in current genome build and once with liftover
    FORMAT_SUMSTATS(
        input_files_ch
            .combine(custom_col_headers)
            .combine(r_lib)
    )

    // emit:
}

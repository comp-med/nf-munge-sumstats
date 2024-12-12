include {
    GET_GENOME_BUILD
} from '../modules/munge_sumstats.nf'

workflow MUNGE_SUMSTATS {

    take: 
    input_files_ch
    custom_col_headers
    r_lib

    main:

    GET_GENOME_BUILD(
        input_files_ch
            .combine(r_lib) 
    ).view()

    // TODO
    // MUNGE_SUMSTATS()

    // emit:
}

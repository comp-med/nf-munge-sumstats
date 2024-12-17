include {
    GET_GENOME_BUILD;
    FORMAT_SUMSTATS;
    SORT_GZIP_INDEX;
} from '../modules/munge_sumstats.nf'

workflow MUNGE_SUMSTATS {

    take: 
    input_files_ch
    custom_col_headers
    r_lib
    bcftools_liftover_bin
    bgzip_bin

    main:

    input_files_ch = GET_GENOME_BUILD (
        input_files_ch.combine(r_lib) 
    ).map {
        tup -> [tup[0], file(tup[1]).text, tup[2]]
    }

    // Format files in their current genome build
    formatted_files_ch = FORMAT_SUMSTATS (
        input_files_ch
            .combine(custom_col_headers)
            .combine(r_lib)
    )

    // Before liftover, files need to be gzipped and indexed
    indexed_files_ch = SORT_GZIP_INDEX (
        formatted_files_ch
            .combine(bcftools_liftover_bin)
            .combine(bgzip_bin)
    ).view()

    //def liftover_files_ch = GET_LIFTOVER_FILES (
        // TODO
    //)

    //def liftover_ch = LIFTOVER_SUMSTATS (
        // TODO
    //)

    // TODO: Join formatted files in both genome builds as output
    // output formatted data directly, then output liftover files

    // emit:
}

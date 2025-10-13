// Check input files for custom headers
// Create the custom header file for mapping column names

include {
    CHECK_INPUT_COL_HEADERS;
    GET_INPUT_COL_HEADERS
} from '../modules/setup_munging.nf'

workflow SETUP_MUNGING {

    take: 
    input_files_ch
    r_lib

    main:

    input_files_ch = input_files_ch
        // Do NOT match VCF(.gz) files
        .filter({tuple -> tuple[1] =~ /^(?!.*\.vcf(?:\.gz)?$).+/})

    input_files_ch = GET_INPUT_COL_HEADERS(
        input_files_ch.combine(r_lib)
    ).map { file -> "$file"
    }.collectFile(
        name: "all_input_files.csv",
        newLine: true
    )

    def check_input_col_headers_ch = CHECK_INPUT_COL_HEADERS (
        input_files_ch,
        r_lib
    )

    emit:
    custom_col_headers = check_input_col_headers_ch
}

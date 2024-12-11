// Check input files for custom headers
// Create the custom header file for mapping column names

include {
    CHECK_INPUT_COL_HEADERS 
} from '../modules/setup_munging.nf'

workflow SETUP_MUNGING {

    take: 
    input_files_ch
    r_lib

    main:
    def input_file_table_ch = input_files_ch
        .map {
            row -> "${row[1]}"
        }
        .collectFile(
            name: "all_input_files",
            newLine: true
        )

    def custom_col_headers = CHECK_INPUT_COL_HEADERS (
        input_file_table_ch,
        r_lib
    )

    custom_col_headers.view()



    // emit:
}

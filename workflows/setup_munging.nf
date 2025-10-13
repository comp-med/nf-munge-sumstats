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
    // TODO: Create a new Process to read column names to then pass into CHECK...
    def input_file_table_ch = input_files_ch
        .map {
            row -> "${row[1]}"
        }
        .collectFile(
            name: "all_input_files",
            newLine: true
        )

    def check_input_col_headers_ch = CHECK_INPUT_COL_HEADERS (
        input_file_table_ch,
        r_lib
    )

    emit:
    custom_col_headers = check_input_col_headers_ch
}

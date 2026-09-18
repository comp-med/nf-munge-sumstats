// Check input files for custom headers
// Create the custom header file for mapping column names

include {
    CHECK_INPUT_COL_HEADERS;
    GET_INPUT_COL_HEADERS;
    DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES
    SETUP_LIFTOVER
    CHECK_LIFTOVER_PLUGIN
} from '../modules/setup_munging.nf'

workflow SETUP_MUNGING {

    take: 
    input_files_ch
    corrected_column_names
    redundant_column_names

    main:

    input_files_ch = input_files_ch
        // Do NOT match VCF(.gz) files
        .filter({tuple -> tuple[1] =~ /^(?!.*\.vcf(?:\.gz)?$).+/})

    input_files_ch = GET_INPUT_COL_HEADERS(
        input_files_ch
    ).map { file -> "$file"
    }.collectFile(
        name: "all_input_files.csv",
        newLine: true
    )

    CHECK_INPUT_COL_HEADERS (
        input_files_ch
          .combine(corrected_column_names)
          .combine(redundant_column_names)
    )

    DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES()

    // Install the liftover plugin in the stock conda env and check it
    SETUP_LIFTOVER()
    CHECK_LIFTOVER_PLUGIN(
        SETUP_LIFTOVER.out
    )

    emit:
    col_headers = CHECK_INPUT_COL_HEADERS.out
    snplocs_lib = DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES.out[0] // Only base dir

}

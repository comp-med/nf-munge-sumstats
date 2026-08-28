// Check input files for custom headers
// Create the custom header file for mapping column names

include {
    CHECK_INPUT_COL_HEADERS;
    GET_INPUT_COL_HEADERS;
    DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES
} from '../modules/setup_munging.nf'

workflow SETUP_MUNGING {

    take: 
    input_files_ch

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
    )

    DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES()

    emit:
    col_headers = CHECK_INPUT_COL_HEADERS.out
    snplocs_lib = DOWNLOAD_BIOCONDUCTOR_DEPENDENCIES.out

}

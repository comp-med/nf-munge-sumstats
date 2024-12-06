// Download Summary statistics data from different sources and combine
// them into a single channel with a phenotype name and a local file link.

include { 
    DOWNLOAD_OTHER_DATA;
    DOWNLOAD_OPEN_GWAS_DATA;
    GWAS_CATALOG_SETUP;
    DOWNLOAD_GWAS_CATALOG_DATA;
    } from '../modules/download_data.nf'

workflow DOWNLOAD_DATA {
    /* 
    * Read content from input file and parse it into separate channels based on
      input source
    * Create separate channels for each kind of data source
    * Also clean input data based on what is needed (ID, link, location)
    */
    take:
    input_table
    r_lib

    main:
    def data_src_ch = Channel.fromPath(
        input_table,
        checkIfExists: true,
        type: 'file'
        ).splitCsv( header: true )
            .map { row -> [
                row.phenotype_id,
                row.data_source, 
                row.data_id, 
                row.data_link, 
                row.data_location
            ] }
            .branch { tup -> 
                gwas_catalog: tup[1] == 'gwas_catalog'
                    return [tup[0], tup[2]]
                open_gwas: tup[1] == 'open_gwas'
                    return [tup[0], tup[2]]
                other: tup[1] == 'other'
                    return [tup[0], tup[3]]
                local: tup[1] == 'local'
                    return [tup[0], file( tup[4] )]
            }

    // Process each channel separately based on sub-channel
    def other_data_ch = DOWNLOAD_OTHER_DATA( data_src_ch.other )
    def open_gwas_data_ch = DOWNLOAD_OPEN_GWAS_DATA( data_src_ch.open_gwas )

    def gwas_cat_setup_ch = GWAS_CATALOG_SETUP( r_lib )
    def gwas_cat_data_ch = DOWNLOAD_GWAS_CATALOG_DATA( 
        data_src_ch.gwas_catalog
            .combine(gwas_cat_setup_ch)
            .combine(r_lib)
    )
    
    // Join channels that contain the sumstat files
    def data_ch = data_src_ch
        .local
        .concat(
            other_data_ch,
            gwas_cat_data_ch,
            open_gwas_data_ch
    )

    emit:
    data = data_ch
}

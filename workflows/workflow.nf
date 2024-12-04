// include { PROCESS_NAME } from '../modules/module_template.nf'

workflow PARSE_INPUT {
  take:
  input_table

  main:


  //  data_id data_link data_location
  Channel.fromPath(
      input_table,
      checkIfExists: true,
      type: 'file')
  .splitCsv( header: true )
  .map {row -> [
    row.phenotype_id,
    row.data_source, 
    row.data_id, 
    row.data_link, 
    row.data_location
  ]}
  .view()

}

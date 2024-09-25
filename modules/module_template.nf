process PROCESS_NAME {
  
  // cache 'lenient'
  // tag "$input_value"
  // label 'process_label'
  // publishDir

  input:
  // tuple val(some_value), path(some_path)
  // path(some_path)
  // each path(some_path)

  output:
  // path "./optional_dir/*file.ext", optional: false, emit: output_name

  script:
  """
  ./script.R \
    $input_parameter \
    ${task.some_parameter} \
    ${params.some_parameter}

  """

}

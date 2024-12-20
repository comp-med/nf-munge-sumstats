#!/usr/bin/env nextflow

// https://sydney-informatics-hub.github.io/template-nf-guide/
nextflow.enable.dsl=2

// HEADER ---------------------------------------------------------------------

log.info """\
===============================================================================
nf-Munge-Sumstats
===============================================================================

Created by the Computational Medicine Group | BIH @ Charité

===============================================================================
Workflow run parameters 
===============================================================================
input       : ${params.input}
outDir      : ${params.outDir}
workDir     : ${workflow.workDir}
===============================================================================

"""

// Help function
def helpMessage() {
  log.info"""
  Usage:  nextflow run main.nf 

  Required Arguments:

  <TODO>

  Optional Arguments:

  --outDir	Specify path to output directory. Default is `output/`
	
""".stripIndent()
}

// WORKFLOWS ------------------------------------------------------------------

include { SETUP_MUNGING } from './workflows/setup_munging.nf'
include { MUNGE_SUMSTATS } from './workflows/munge_sumstats.nf'

// WORKFLOW -------------------------------------------------------------------

workflow {

  def input_dir = file("$params.input")

  // Where to find all R packages
  def r_lib    = Channel.fromPath(params.local_r_library)

  // Where to find additional binaries // TODO: create environments
  def bcftools_liftover_bin = Channel.fromPath(params.bcftools_liftover_bin)
  def bgzip_bin = Channel.fromPath(params.bgzip_bin)
  
  // Download raw summary statistics from various sources
  def input_files_ch  = Channel
      .fromPath(
          "$input_dir/**/raw_sumstat_file.*",
          followLinks: true,
          checkIfExists: true)
       .map { 
      path -> [path.getParent().getName(), file(path)]
  }

  // Prepare additional input for liftover function
  SETUP_MUNGING (input_files_ch, r_lib)
  def custom_col_headers = SETUP_MUNGING.out.custom_col_headers

  // Main workflow: format and liftover summary statistics
  def munged_sumstats_ch = MUNGE_SUMSTATS(
      input_files_ch,
      custom_col_headers,
      r_lib,
      bcftools_liftover_bin,
      bgzip_bin
  )
}

// SUMMARY --------------------------------------------------------------------

workflow.onComplete {
    def summary = """\
===============================================================================
Workflow execution summary
===============================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.outDir}

===============================================================================
"""
    println summary
}


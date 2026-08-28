#!/usr/bin/env nextflow

// https://sydney-informatics-hub.github.io/template-nf-guide/
nextflow.enable.dsl=2

// HEADER ---------------------------------------------------------------------

def startupMsg() {
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

  """.stripIndent()
}

def completionMsg() {
  log.info """\
  ===============================================================================
  Workflow execution summary
  ===============================================================================

  Duration    : ${workflow.duration}
  Success     : ${workflow.success}
  workDir     : ${workflow.workDir}
  Exit status : ${workflow.exitStatus}
  outDir      : ${params.outDir}

  ===============================================================================
  """.stripIndent()
}

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

// ENTRY WORKFLOW -------------------------------------------------------------

workflow {

  startupMsg()
  if ( false ) {
    helpMessage()
    exit 1
  }

  // Create Variables from parameters // TODO: Have all parameters here!
  def input_dir = file("$params.input")

  // Where to find additional binaries // TODO: create environments
  def bcftools_liftover_bin = channel.fromPath(params.bcftools_liftover_bin)
  
  // Download raw summary statistics from various sources
  def input_files_ch  = channel
      .fromPath(
          "$input_dir/**/raw_sumstat_file.*",
          followLinks: true,
          checkIfExists: true)
       .map { 
      path -> [path.getParent().getName(), file(path)]
  }

  // Prepare additional input for liftover function
  SETUP_MUNGING (input_files_ch)

  def custom_col_headers = SETUP_MUNGING.out.col_headers
  def snplocs_lib = SETUP_MUNGING.out.snplocs_lib

  // Main workflow: format and liftover summary statistics
  MUNGE_SUMSTATS(
    input_files_ch,
    custom_col_headers,
    snplocs_lib,
    bcftools_liftover_bin
  )

  workflow.onComplete {
      completionMsg()
  }

}


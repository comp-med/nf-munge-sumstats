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

include { DOWNLOAD_DATA } from './workflows/download_data.nf'
include { SETUP_MUNGING } from './workflows/setup_munging.nf'
include { MUNGE_SUMSTATS } from './workflows/munge_sumstats.nf'

// WORKFLOW -------------------------------------------------------------------

workflow {

  // Main file containing Data IDs and necessary meta data
  def input_table = file(params.input_table)

  // Where to find all R packages
  def r_lib    = Channel.fromPath(params.local_r_library)

  // Where to find additional binaries // TODO: create environments
  def lftp_bin = Channel.fromPath(params.lftp_bin)
  def bcftools_liftover_bin = Channel.fromPath(params.bcftools_liftover_bin)
  def bgzip_bin = Channel.fromPath(params.bgzip_bin)
  
  // Download raw summary statistics from various sources
  DOWNLOAD_DATA (input_table, r_lib, lftp_bin)
  def input_files_ch = DOWNLOAD_DATA.out.data

  // Prepare additional input for liftover function
  SETUP_MUNGING (input_files_ch, r_lib)
  def custom_col_headers = SETUP_MUNGING.out.custom_col_headers

  // Main workflow: format and liftover summary statistics
  MUNGE_SUMSTATS(
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


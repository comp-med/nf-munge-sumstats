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

// PLAN -----------------------------------------------------------------------

/*

* Parse input table - what do I need as input information?
* Download data
* Run it through munge sumstats
* Liftover

* Automate GWAS Catalog and openGWAS downloads
* For everything else, have download link ready
* Distinguish between VCF and regular table
* Save VCF as tsv.gz in separate step

Steps:

1. Compile Data
 - Define mandatory input columns: phenotype_name, download_source [gwas_catalog, open_gwas, other], download_link, local File [T, F], file_location
2. Download Data
 - Split by download type, skip locally available data
3. Prepare harmonization 
 - create custom `sumstatsColHeaders` (How to output unmatched colnames?)
 - Get genome build for each file irrespective of annotation, which can be wrong
 - 
3. Harmonize Data (Both Genome Builds)

Support: GBMI, Zenodo

*/

// WORKFLOWS ------------------------------------------------------------------

include { DOWNLOAD_DATA } from './workflows/download_data.nf'
include { SETUP_MUNGING } from './workflows/setup_munging.nf'
include { MUNGE_SUMSTATS } from './workflows/munge_sumstats.nf'

// WORKFLOW -------------------------------------------------------------------

workflow {

  // Main file containing Data IDs and necessary meta data
  def input_table = file(params.input_table)

  // Where to find all R packages
  def r_lib = Channel.fromPath(params.local_r_library)
  def lftp_bin = Channel.fromPath(params.lftp_bin)
  
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
    r_lib
  )
  // Check integrity of files somehow?
  // Create colheaders table
  // DEBUG output
  // input_files_ch.view()

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


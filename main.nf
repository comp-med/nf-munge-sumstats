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
3. Harmonize Data (Both Genome Builds)

Support: GBMI, Zenodo

*/

// MODULES --------------------------------------------------------------------

include { DOWNLOAD_DATA } from './workflows/download_data.nf'

// WORKFLOW -------------------------------------------------------------------

workflow {

  // Main file containing Data IDs and necessary meta data
  def input_table = file(params.input_table)

  // Where to find all R packages
  def r_lib = Channel.fromPath(params.local_r_library)
  
  // Download raw summary statistics from various sources
  DOWNLOAD_DATA (input_table, r_lib)

  // DEBUG
  DOWNLOAD_DATA.out.data.view()

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


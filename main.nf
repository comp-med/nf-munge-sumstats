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

*/

// MODULES --------------------------------------------------------------------

include { WORKFLOW } from './workflows/workflow.nf'

// WORKFLOW -------------------------------------------------------------------

workflow {
  
  WORKFLOW ()

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


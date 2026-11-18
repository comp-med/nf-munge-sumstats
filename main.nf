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

  --input         Directory containing the raw summary statistic files

  Optional Arguments:

  --outDir	  Specify path to output directory. Default is `output/`
  --proxy_server  When running pipeline in a network with a proxy, specify here
	
""".stripIndent()
}

// WORKFLOWS ------------------------------------------------------------------

include { SETUP_MUNGING } from './workflows/setup_munging.nf'
include { MUNGE_SUMSTATS } from './workflows/munge_sumstats.nf'

// ENTRY WORKFLOW -------------------------------------------------------------

workflow {

  main:

  if (params.verbose == true ) startupMsg()
  if ( params.input == null ) {
    helpMessage()
    exit 1
  }

  def input_dir = file("$params.input")

  def corrected_column_names = channel.fromPath(params.corrected_column_names)
  def redundant_column_names = channel.fromPath(params.redundant_column_names)
  def log_pvalue_column_names = channel.fromPath(params.log_pvalue_column_names)

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
  SETUP_MUNGING (
    input_files_ch,
    corrected_column_names,
    redundant_column_names,
    log_pvalue_column_names
  )

  def custom_col_headers = SETUP_MUNGING.out.col_headers
  def snplocs_lib = SETUP_MUNGING.out.snplocs_lib

  // Main workflow: format and liftover summary statistics
  MUNGE_SUMSTATS(
    input_files_ch,
    custom_col_headers,
    snplocs_lib,
    log_pvalue_column_names
  )

  if (params.verbose == true ) {
    workflow.onComplete {
        completionMsg()
    }
  }
  
  publish:
    parquet = MUNGE_SUMSTATS.out.parquet
    check_rejected = MUNGE_SUMSTATS.out.check_rejected
    vcf = MUNGE_SUMSTATS.out.vcf
    format_log = MUNGE_SUMSTATS.out.format_log
    liftover_log = MUNGE_SUMSTATS.out.liftover_log
    liftover_rejected = MUNGE_SUMSTATS.out.liftover_rejected

}

// OUTPUTS -------------------------------------------------------------

output {
      parquet {
        path { pheno, build, file -> 
        [ (file): "formatted/${pheno}/${build}/formatted_sumstats.parquet" ] 
        }
        enabled params.publish_parquet
        mode 'copy'
    }
    vcf {
        path { pheno, build, file, index -> [
           (file): "formatted/${pheno}/${build}/",
           (index): "formatted/${pheno}/${build}/" 
           ] }
        enabled params.publish_vcf
        mode 'copy'
    }
    format_log {
        path { pheno, file -> [
          (file): "formatted/${pheno}/logs/formatting_log.txt"
        ] }
        mode 'copy'
    }
    liftover_log {
        path { pheno, file -> [
          (file): "formatted/${pheno}/logs/",
        ] }
        mode 'copy'
    }
    liftover_rejected {
        path { pheno, file -> [
          (file): "formatted/${pheno}/rejected/liftover_rejected.vcf.bgz"
        ] }
        mode 'copy'
    }
    check_rejected {
        path { pheno, file -> [
          (file): "formatted/${pheno}/rejected/check_rejected.parquet"
        ] }
        mode 'copy'
    }

}
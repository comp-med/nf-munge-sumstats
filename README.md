# nf-munge-sumstats

## Introduction

This is a Nextflow Pipeline to streamline formatting and liftover of GWAS
summary statistics into a harmonized format
([GWAS-VCF](https://github.com/MRCIEU/gwas-vcf-specification)). This pipeline
is intended to be used with the output of the companion pipeline
[`nf-download-sumstats`](https://github.com/cfbeuchel/nf-download-sumstats).

The pipeline covers the following steps:

1. Format summary statistics according to GWAS-VCF standard using the
   [`MungeSumstats`](https://academic.oup.com/bioinformatics/article/37/23/4593/6380562)
   R-package
2. Based on raw data reference genome alignment, liftover formatted file to the
   other reference genome (either GRCh37 or GRCh38) using
   [`bcftools\liftover`](https://academic.oup.com/bioinformatics/article/40/2/btae038/7585532)
3. Save files as indexed VCF files as well as parquet tables


## Requirements

### Input Data: Raw Data Directory 

This pipeline requires a structured directory containing GWAS summary
statistics files that was generated by the Nextflow pipeline
[`nf-download-sumstats`](https://github.com/cfbeuchel/nf-download-sumstats).
See the repository for the required data structure. In short, for each summary
statistics file, the phenotype name will be parsed from the parent directory
and the file containing the data must be named `raw_sumstat_file`, with a file
extension, e.g. `.gz`, `.vcf`. Based on the default output of the downloading
pipeline, the `input` parameter should point to the `raw/` directory.

### Software

Currently, the pipeline does not automatically set up all required software and
R packages so make sure you have the following software dependencies set up and
configured in `nextflow.config`.

* External binaries
    * [`bgzip`](https://www.htslib.org/doc/bgzip.html)
    * [`bcftools`](https://samtools.github.io/bcftools/bcftools.html) with the [`liftover`](https://github.com/freeseek/score/?tab=readme-ov-file#liftover-vcfs) plugin installed
* R packages: 
    * From GitHub: [MungeSumstats](https://github.com/Al-Murphy/MungeSumstats)
    * From CRAN: `data.table`, `fs`, `arrow`
    * From Bioconductor: `GenomicFiles`, `GenomicRanges`, `VariantAnnotation`,
      `GenomeInfoDb`, `SNPlocs.Hsapiens.dbSNP155.GRCh37`, 
      `SNPlocs.Hsapiens.dbSNP155.GRCh38`

## Gettings Started

Please check out the
[`nf-download-sumstats`](https://github.com/cfbeuchel/nf-download-sumstats)
pipeline on how to download raw summary statistics and then use the output
produced by that pipeline as the `input` for this pipelin.

Make sure to add all the required parameters and path to the input directory in
`nextflow.config`.  Afterwards, you can run the pipeline.

First, see if the the pipeline executes by starting a dry-run.

```bash
# When running the pipeline locally, e.g. on your laptop
nextflow run main.nf -stub-run 

# When on a HPC with SLURM, use the `cluster` profile
nextflow run main.nf -stub-run -profile cluster
```

If that finished successfully, run the actual pipeline.

```bash
# Again, when on an HPC with SLURM, add the `-profile cluster` flag
nextflow run main.nf -profile cluster
```

### Output

By default, all output will be saved in the `./output` directory in the
pipeline's base directory. Formatted files will be copied into the
`./formatted` sub-directory. Files for each phenotype will be saved in
subdirectories named after each line in `phenotype_id` of `params.input_table`.

```bash
output/
├── formatted
│   ├── phentype_1
│   │   ├── grch37
│   │   │   ├── formatted_sumstats_grch37.parquet
│   │   │   ├── formatted_sumstats_grch37.vcf.gz
│   │   │   └── formatted_sumstats_grch37.vcf.gz.tbi
│   │   └── grch38
│   │       ├── formatted_sumstats_grch38.parquet
│   │       ├── formatted_sumstats_grch38.vcf.gz
│   │       └── formatted_sumstats_grch38.vcf.gz.tbi
│   ├── phenotype_2
│   │   ├── grch37
│   │   │   ├── formatted_sumstats_grch37.parquet
│   │   │   ├── formatted_sumstats_grch37.vcf.gz
│   │   │   └── formatted_sumstats_grch37.vcf.gz.tbi
│   │   └── grch38
│   │       ├── formatted_sumstats_grch38.parquet
│   │       ├── formatted_sumstats_grch38.vcf.gz
│   │       └── formatted_sumstats_grch38.vcf.gz.tbi
[...]
```

## Expected Problems

The pipeline will fail when an input file contains column headers that are not
covered by the mapping table provided by `MungeSumstats` and augmented in the
`CHECK_INPUT_COL_HEADERS` process. All column names need to be mapped in this
step, either as those that can be corrected or as superfluous. You might need
to manually adjust the code in this step to make it run successfully.

##  TODO

* Make more setting for `MungeSumstats` user-facing
* Externalize the currently internal proxy settings
* Look for `TODO`s in the code!
* Make propagating meta data more ergonomic
* Add support for downloading files from Zenodo (& GBMI?)
* More fine-grained resource allocation
* Add nf-test
* Add support for more reporting
* Add help text to command
* Create nice output table with with phenotype & path

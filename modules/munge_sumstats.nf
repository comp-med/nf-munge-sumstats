// To know whether I want to liftover data or not I need to infer genome build
process GET_GENOME_BUILD {
    
    cache 'lenient'
    tag "$phenotype_name"
    label 'rProcess'
    conda 'conda_envs/r.yml'

    input:
    tuple val(phenotype_name), 
        path(raw_sumstat_file), 
        path(custom_col_headers),
        path(snplocs_lib)

    output:
    tuple val(phenotype_name), path("genome_build"), path(raw_sumstat_file)

    script:
    """
    get_genome_build.R \
        --pheno $phenotype_name \
        --file $raw_sumstat_file \
        --col_headers $custom_col_headers \
        --lib $snplocs_lib \
        --cores $task.cpus \
    """

    stub:
    """
    echo "grch37" > genome_build
    """

}

// Format summary statistics without lifting them over
process FORMAT_SUMSTATS {
    
    cache 'lenient'
    tag "$phenotype_name, $genome_build"
    label 'rProcess'
    conda 'conda_envs/r.yml'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path(raw_sumstat_file), 
        path(custom_col_headers),
        path(snplocs_lib),
	    path(log_pvalue_column_names_file)

    output:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path("formatted_sumstats.vcf"), emit: files
    tuple val(phenotype_name), path("formatted_sumstats_log_msg.txt"), emit: log

    script:
    """
    format_sumstats.R --pheno $phenotype_name \
        --build $genome_build \
        --file $raw_sumstat_file \
        --col_headers $custom_col_headers \
        --lib $snplocs_lib \
        --logp_names $log_pvalue_column_names_file
    """

    stub:
    """
    touch "formatted_sumstats.vcf"
    touch "formatted_sumstats_log_msg.txt"
    """

}

// Get chain files and reference sequences necessary for liftover
process GET_LIFTOVER_FILES {
    
    cache 'lenient'
    tag 'single_execution'
    conda 'conda_envs/samtools.yml'

    output:
    tuple path('hg19.fa'),
        path('hg38.fa'),
        path('hg19ToHg38.over.chain.gz'),
        path('hg38ToHg19.over.chain.gz')

    script:
    """
   
    # Chain Files
    CHAIN1='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz'
    CHAIN2='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz'
    wget \$CHAIN1
    wget \$CHAIN2

    # Reference Assemblies
    REF1='https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz'
    REF2='https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz'
    wget \$REF1
    wget \$REF2

    # Unzip the reference sequences and remove the the 
    bgzip -d  \$(echo "\${REF1##*/}")
    bgzip -d  \$(echo "\${REF2##*/}")
    """

    stub:
    """
    touch hg19ToHg38.over.chain.gz hg38ToHg19.over.chain.gz hg19.fa hg38.fa
    """

}

// When file is already in output format, do not perform liftover
process NO_LIFTOVER_SUMSTATS {
    cache 'lenient'
    tag "$phenotype_name:$genome_build"
    conda 'conda_envs/bcftools_liftover.yml'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path(formatted_sumstats),
        path(hg19_reference),
        path(hg38_reference),
        path(_hg19_to_38_chain_file),
        path(_hg38_to_19_chain_file)

    output:
    tuple val(phenotype_name),
        val(genome_build),
        path("formatted_sumstats_${genome_build}.vcf.gz"),
        emit: files
    tuple val(phenotype_name), 
        path("no_liftover_log.txt"),
        emit: log

    script:
    if (!(genome_build in ['grch37', 'grch38']))
    error "Unknown genome_build: ${genome_build}"
    def forward = genome_build == 'grch37'
    def src_ref = forward ? hg19_reference : hg38_reference
    """
    no_liftover.sh \
        --input   ${formatted_sumstats} \
        --from    ${genome_build} \
        --src-ref ${src_ref}
    """
    stub:
    """
    touch "formatted_sumstats_{$genome_build}.vcf.gz"
    touch "no_liftover.log"
    """
}

// Liftover summary statistics from given assembly to the other one
process LIFTOVER_SUMSTATS {

    cache 'lenient'
    tag "$phenotype_name:$genome_build->$other_genome_build"
    conda 'conda_envs/bcftools_liftover.yml'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path(formatted_sumstats),
        path(hg19_reference),
        path(hg38_reference),
        path(hg19_to_38_chain_file),
        path(hg38_to_19_chain_file)
    
    output:
    tuple val(phenotype_name),
        val(genome_build),
        val(other_genome_build),
        path("formatted_sumstats_grch37.vcf.gz"),
        path("formatted_sumstats_grch38.vcf.gz"), 
        emit: files
    tuple val(phenotype_name), 
        path("liftover_log.txt"),
        emit: log
    tuple val(phenotype_name), 
        path("reject.vcf.bgz"),
        emit: reject

    script:
    if (!(genome_build in ['grch37', 'grch38']))
    error "Unknown genome_build: ${genome_build}"
    def forward = genome_build == 'grch37'
    def chain   = forward ? hg19_to_38_chain_file : hg38_to_19_chain_file
    def src_ref = forward ? hg19_reference        : hg38_reference
    def tgt_ref = forward ? hg38_reference        : hg19_reference
    """
    liftover.sh \
        --input   ${formatted_sumstats} \
        --from    ${genome_build} \
        --to      ${other_genome_build} \
        --src-ref ${src_ref} \
        --tgt-ref ${tgt_ref} \
        --chain   ${chain}
    """
    
    stub:
    """
    touch "formatted_sumstats_grch37.vcf.gz"
    touch "formatted_sumstats_grch38.vcf.gz"
    touch "liftover_log.txt"
    touch "reject.vcf.bgz"
    """

}

// To make downstream processing easier, save as tabular parquet
process CREATE_PARQUET {

    cache 'lenient'
    tag "$phenotype_name"
    label 'rProcess'
    conda 'conda_envs/r.yml'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        path(formatted_sumstats)

    output:
    tuple val(phenotype_name),
        val(genome_build),
        path("formatted_sumstats.parquet")


    script:
    """
    create_parquet_from_vcf.R $formatted_sumstats
    """

    stub:
    """
    touch formatted_sumstats.parquet
    """

}

process CHECK_PARQUET {

    cache 'lenient'
    tag "$phenotype_name"
    label 'rProcess'
    conda 'conda_envs/r.yml'

    input:
    tuple val(phenotype_name),
        val(genome_build),
        path(formatted_sumstats),
        path(snplocs_lib)


    output:
    tuple val(phenotype_name), 
        val(genome_build), 
        path("checked_sumstats.parquet"),
        emit: parquet
    tuple val(phenotype_name), path("rejected_sumstats.parquet"), 
        optional: true,
        emit: rejected

    script:
    """
    check_parquet.R $formatted_sumstats $snplocs_lib
    """

    stub:
    """
    touch checked_sumstats.parquet
    touch rejected_sumstats.parquet

    """
}

process SAVE_VCF {
    cache 'lenient'
    tag "$phenotype_name"
    label 'rProcess'
    conda 'conda_envs/r.yml'
    
    input:
    tuple val(phenotype_name),
        val(genome_build),
        path(formatted_sumstats),
        path(snplocs_lib)

    output:
    tuple val(phenotype_name), 
        val(genome_build), 
        path("formatted_sumstats.vcf.bgz"),
        path("formatted_sumstats.vcf.bgz.tbi")


    script:
    """
    write_vcf.R $formatted_sumstats $snplocs_lib
    """

    stub:
    """
    touch formatted_sumstats.vcf.bgz
    touch formatted_sumstats.vcf.bgz.tbi
    """


}
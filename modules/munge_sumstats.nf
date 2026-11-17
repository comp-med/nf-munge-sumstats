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
        path(munging_log)
    output:
    tuple val(phenotype_name),
        val(genome_build),
        path("formatted_sumstats.vcf.gz")

    script:
    """
    INPUT_VCF="$formatted_sumstats"
    OUTPUT_VCF="\${INPUT_VCF}.gz"

    # Small pvalues will underflow to 0 since float are 32bit
    # Saving them as String retains the value
    bcftools view -h "\$INPUT_VCF" | \
        sed -E '/^##FORMAT=<ID=P,/ s/Type=Float/Type=String/' > \
        new_header.txt

    # Re-header, sort, expand, index and save
    bcftools reheader -h new_header.txt "\$INPUT_VCF" | \
        bcftools  sort -Ou | \
        bcftools  norm --no-version -Ou -m- | \
        bcftools annotate -x INFO -Oz -o \$OUTPUT_VCF
    """
    stub:
    """
    touch "formatted_sumstats.vcf.gz"
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
        path("reject.vcf.bgz"),
        emit: log

    script:
    """
    # INPUT FILES ----
    HG19_REF="$hg19_reference"
    HG38_REF="$hg38_reference"
    HG19_TO_HG38_CHAIN="$hg19_to_38_chain_file"
    HG38_TO_HG19_CHAIN="$hg38_to_19_chain_file"
    FROM_GENOME_BUILD="$genome_build"
    TO_GENOME_BUILD="$other_genome_build"
    INPUT_VCF="$formatted_sumstats"
    OUTPUT_VCF="formatted_sumstats_\${TO_GENOME_BUILD}.vcf.gz"

    # LOG FILE ----
    touch liftover_log.txt
    LOG='liftover_log.txt'

    echo "Liftover direction: \${FROM_GENOME_BUILD} -> \${TO_GENOME_BUILD}" >> \$LOG 2>&1

    # SORT & BGZIP INPUT ---- 

    # Small pvalues will underflow to 0 since float are 32bit
    # Saving them as String retains the value
    bcftools view -h "\$INPUT_VCF" | \
        sed -E '/^##FORMAT=<ID=P,/ s/Type=Float/Type=String/' > \
        new_header.txt

    # Re-header, sort, expand, index and save
    bcftools reheader -h new_header.txt "\$INPUT_VCF" | \
        bcftools sort -o \$INPUT_VCF >> \$LOG 2>&1

    # BGZip the file
    bgzip \$INPUT_VCF
    INPUT_VCF="\${INPUT_VCF}.gz"
    INPUT_VCF_NAMED="formatted_sumstats_\${FROM_GENOME_BUILD}.vcf.gz"
    mv \$INPUT_VCF \$INPUT_VCF_NAMED
    INPUT_VCF="\${INPUT_VCF_NAMED}"

    # LIFTOVER DIRECTION ----

    if [ "\$FROM_GENOME_BUILD" == "grch37" ]; then
        CHAIN=\$HG19_TO_HG38_CHAIN
        SOURCE_REF=\$HG19_REF
        TARGET_REF=\$HG38_REF
    else
        CHAIN=\$HG38_TO_HG19_CHAIN
        SOURCE_REF=\$HG38_REF
        TARGET_REF=\$HG19_REF    
    fi

    # COLLAPSE VCF ----

    # Create a collapsed VCF and check REF/ALT alignment in the process
    bcftools norm --no-version -Ou -m+ \$INPUT_VCF \
        --check-ref ws \
        -f \$SOURCE_REF \
        -o input_COLLAPSED.vcf.gz \
        >> \$LOG 2>&1

    # FIX VCF REF ----

    # Check Allele mismatches etc
    bcftools +fixref input_COLLAPSED.vcf.gz -- -f \$SOURCE_REF \
        >> \$LOG 2>&1

    # LIFTOVER ----

    # Liftover
    bcftools +liftover --no-version input_COLLAPSED.vcf.gz -Ou -o output_COLLAPSED.vcf.gz -- \
        -s \$SOURCE_REF \
        -f \$TARGET_REF \
        -c \$CHAIN  \
        --reject reject.vcf.bgz \
        --reject-type z \
        --write-src \
        >> \$LOG 2>&1

    # CHECKS ----

    # Check the file using bcftools
    bcftools +af-dist output_COLLAPSED.vcf.gz \
        >> \$LOG 2>&1

    bcftools +fixref output_COLLAPSED.vcf.gz -- -f \$TARGET_REF \
        >> \$LOG 2>&1

    bcftools stats output_COLLAPSED.vcf.gz \
        >> \$LOG 2>&1

    # Sort, expand, index and save
    bcftools  sort -Oz output_COLLAPSED.vcf.gz | \
        bcftools  norm --no-version -Oz -m- -o output.vcf.gz \
        >> \$LOG 2>&1

    # INFO Column of output file contains content that his difficult
    # to parse downstream. I'll remove it
    bcftools annotate -x INFO output.vcf.gz -Oz -o \$OUTPUT_VCF \
        >> \$LOG 2>&1

    # This is because bcftools does not modify in-place, apparently.
    # Clean up
    rm output.vcf.gz
    rm output_COLLAPSED.vcf.gz
    rm input_COLLAPSED.vcf.gz
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
        path("formatted_sumstats.parquet"),
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
    touch formatted_sumstats.parquet
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
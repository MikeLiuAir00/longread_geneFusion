/*This is the Nextflow script for gene fusion detection pipeline with LongGF*/

param.reads = "$projectDir/data/*.{fastq, fq}"
param.ref = "$projectDir/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa "
param.gtf = "$projectDir/ref/Homo_sapiens.GRCh38.109.gtf"
param.result_path= "$projectDir/result"
param.log_path = "$projectDir/log"

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    referecne    : ${params.ref}
    GTF          : ${params.gtf}
    reads        : ${params.reads}
    outdir       : ${params.result_path}
    log          :${param.log_path}
    """
    .stripIndent()


process INDEX_REF {
    // samtools
    container 'quay.io/biocontainers/samtools:1.18--hd87286a_2'
    tag "Indexing $reference"

    input:
    path reference

    output:
    'samtools_faidx'

    script:
    """
    samtools faidx $reference
    """
}

process RAWSTAT {
    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'

    input:
    path read
    path outdir

    output:

}

process PORECHOP {
    container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'
}

process QC {
    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'
}

process MAPPING {
    container 'quay.io/biocontainers/minimap2:2.24--he4a0461_2'
}

process SORTING {
    container 'uay.io/biocontainers/samtools:1.13--hd87286a_2'
}

process FUSION_DETECT {
    container 'quay.io/biocontainers/longgf:0.1.2--h84372a0_6'
}

workflow {
}

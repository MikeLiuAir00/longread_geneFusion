params.reads = "$projectDir/data/*.{fq, fastq}"
params.ref = "$projectDir/ref/hg38.fa"
params.gtf= "$projectDir/ref/hg38.gtf"
params.outdir = "$projectDir/result"
params.out = "$projectDir"


log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    log          : ${params.out}
    """
    .stripIndent()

process REF_INDEX {
    conda "nanoporetech::samtools=1.18"
    tag "indexing reference genome $ref"

    input:
    tuple val(meta), path(ref)

    output:
    tuple val(meta), path('*.fai')

    script:
    """
    samtools faidx $ref
    """

}

process QCSTAT {
    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'
    conda 'bioconda::nanostat=1.6.0'
    tag "stat on $read"
    publishDir "$params.outdir/nanostat", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path ("*.xls")

    // no output dir is specified since each Nextflow process has its own dir(env)
    script:
    """
    NanoStat -t 6 \
        --fastq $read \
        -n ${meta.sample_id}.rawstat.xls
    """
}

process PORECHOP {
    container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'
    conda 'nanoporetech::porechop=0.2.4'
    tag "start porechop trimming adaptor"
    publishDir "$params.outdir/porechop", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path ('*.fastq')

    script:
    """
    porechop -t 6 \
        -i $read \
        -o ${meta.sample_id}_deadapter.fastq \
        --format fastq
    """
}

process NANOFILT {
    // container
    conda 'bioconda::nanofilt=2.8.0'
    tag "Start QC filtering on $read"
    publishDir "$params.outdir/nanofilt", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path('*.gz')

    """
    NanoFilt -q 7 -l 100 \
        $read |
    gzip > ${meta.sample_id}.filtered.fastq.gz
    """
}


process POST_QCSTAT {
    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'
    conda 'bioconda::nanostat=1.6.0'
    tag "stat on $read"
    publishDir "$params.outdir/nanostat", mode: 'copy'

    input:
    tuple val(meta), path(read)

    output:
    tuple val(meta), path ("*.xls")

    // no output dir is specified since each Nextflow process has its own dir(env)
    script:
    """
    NanoStat -t 6 \
        --fastq $read \
        -n ${meta.sample_id}.filtered.stat.xls
    """
}

// process MINIMAP2 {
//     conda 'bioconda::minimap2=2.24'
//     tag 'mapping $read to ref $ref'
// }

workflow {
    // init input channel
    // fastq file input
    Channel.fromPath(params.reads)
    | map {file ->
        sample_id = file.name.tokenize('.')[0]
        meta = [sample_id: sample_id]
        [meta, file]
    }
    | set { fastq_read }

    // reference input
    Channel.fromPath(params.ref)
    | map{ file->
        ref_id = file.name.tokenize('.')[0]
        meta = [ref_id: ref_id]
        [meta, file]
    }
    | set { reference }


    // Workflow start
    REF_INDEX(reference)
    QCSTAT(fastq_read)
    PORECHOP(fastq_read)
    NANOFILT(PORECHOP.out)
    POST_QCSTAT(NANOFILT.out)
}

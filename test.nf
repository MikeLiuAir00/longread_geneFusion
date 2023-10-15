params.reads = "$projectDir/data/*.{fq, fastq}"
params.outdir = "$projectDir/result"
params.out = "$projectDir"
fasta_read = Channel.fromPath(params.reads)
                    .map{tuple(it.name.tokenize('.')[0], it)}

log.info """\
    R N A S E Q - N F   P I P E L I N E
    ===================================
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    log          : ${params.out}
    """
    .stripIndent()


process RAWSTAT {
    container 'quay.io/biocontainers/nanostat:1.6.0--pyhdfd78af_0'
    conda 'bioconda::nanostat=1.6.0'
    tag "stat on $read"
    publishDir "$params.outdir/nanostat", mode: 'copy'

    input:
    tuple val(sample_id), path(read)

    output:
    path "*.xls"

    // no output dir is specified since each Nextflow process has its own dir(env)
    script:
    """
    NanoStat -t 6 \
        --fastq $read \
        -n ${sample_id}.rawstat.xls
    """
}

process PORECHOP {
    container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'
    conda 'nanoporetech::porechop=0.2.4'
    tag "start porechop trimming adaptor"

    input:
    tuple val(sample_id), path(read)

    output:
    path '*.fastq'

    script:
    """
    porechop -t 6 \
        -i $read \
        -o ${sample_id}_deadapter.fastq \
        --format fastq
    """
}

workflow {
    nano_out = RAWSTAT(fasta_read)
    porechop_out = PORECHOP(fasta_read)
    nano_out.view()
    porechop_out.view()
}

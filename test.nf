params.reads = "$projectDir/data/*.{fq, fastq}"
params.outdir = "$projectDir/result"
params.out = "$projectDir"
fasta_read = Channel.fromPath(params.reads)

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
    path read

    output:
    path "*.xls"

    // no output dir is specified since each Nextflow process has its own dir(env)
    script:
    """
    NanoStat -t 6 \
        --fastq $read \
        -n testout.xls
    """
}


workflow {
    nano_out = RAWSTAT(fasta_read)
    nano_out.view()
}

// PARAMETERS
params.CONTAINER = "theoaphidian/wgs-illumina"
params.OUTPUT = "results"

params.MINLEN = 50
params.TRIMLEFT = 0 
params.TRIMRIGHT = 0
params.TRUNCLEN = 0
params.MAXN = 2

process FASTP {
    tag "${pair_id}"    
    publishDir "${params.OUTPUT}/fastp", mode: 'copy'
    
    container params.CONTAINER

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("filtered_${pair_id}*.fastq.gz"), emit: filteredReads
    path("${pair_id}_fastp.json"), emit: fastp
    path("versions.yml"), emit: versions

    script:
    def single = reads instanceof Path

    def input = !single ? "-i '${reads[0]}' -I '${reads[1]}'" : "-i '${reads}'"
    def output = !single ? "-o 'filtered_${pair_id}_fwd.fastq.gz' -O 'filtered_${pair_id}_rev.fastq.gz'" : "-o 'filtered_${pair_id}.fastq.gz'"

    """
    fastp ${input} ${output} --json ${pair_id}_fastp.json \\
     --length_required ${params.minLen} --trim_front1 ${params.trimLeft} \\
     --trim_tail1 ${params.trimRight} --max_len1 ${params.truncLen} \\
     --n_base_limit ${params.maxN} 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastp: \$(fastp -v 2>&1 > /dev/null | sed 's/fastp //';)
    END_VERSIONS    
    """
}

process MULTIQC {
    publishDir "${params.OUTPUT}", mode: 'move', pattern: '*.html'

    container params.CONTAINER
    
    input:
    path('fastp/*')

    output:
    path("multiqc_report.html")
    path("versions.yml"), emit: versions

    script:
    """
    multiqc -m fastp .

        cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //';)
    END_VERSIONS   
    """
}

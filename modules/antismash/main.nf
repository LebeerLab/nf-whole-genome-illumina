params.annotation_folders = "${launchDir}/results/run01/*/annotation"
params.RUNNAME="run01"
params.OUTDIR="results"

process ANTISMASH {
    container null
    tag "${pair_id}"

    publishDir "${params.OUTDIR}/${params.RUNNAME}/${pair_id}", mode: 'move', pattern: "antismash/*"

    input:
    tuple val(pair_id), path(annotation) 

    output:
    tuple val(pair_id), path("antismash/*"), emit: results
    path("versions.yml"), emit: versions
    script:
    """
    gunzip -c $annotation/*.gbff.gz > antismash.gbff
    run_antismash antismash.gbff . -c ${task.cpus} \
      --genefinding-tool none

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(docker run antismash/standalone . --help | grep -oP '\\d{1,2}.\\d{1,2}.\\d{1,2}';)
    END_VERSIONS    
    """

}

workflow {

    annotations = channel.fromPath(params.annotation_folders, type:'dir')
    annotations_id = annotations.map { tuple(file(it).getParent().getName(), it) }
    ANTISMASH(annotations_id)

}

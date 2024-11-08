process SPADES_ASSEMBLY {
    errorStrategy 'ignore'
    container "staphb/shovill:latest"

    tag "${pair_id}" 
    label 'big_mem'

    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("assembly"), emit: assembly
    path("versions.yml"), emit: versions

    script:
    def single = reads instanceof Path

    def input = !single ? "--R1 '${reads[0]}' --R2 '${reads[1]}'" : "--R1 '${reads}'"
    """
    shovill --outdir assembly ${input} --ram ${task.memory} --cpus ${task.cpus} \
        --depth ${params.depth} \
        --minlen ${params.minLenContig} \
        --mincov ${params.minCov} \
	--namefmt "${pair_id}-contig%05d"
    mv assembly/contigs.fa "assembly/${pair_id}_contigs.fna"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        shovill: \$(shovill -v | grep -oP "\\d{1,2}.\\d{1,2}.\\d{1,2}" ;)
    END_VERSIONS    
    """
}

process FLYE_ASSEMBLY {
errorStrategy 'ignore'
    container "staphb/flye:latest"

    tag "${pair_id}" 
    label 'big_mem'

    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("assembly"), emit: assembly
    path("versions.yml"), emit: versions

    script:
    """
    flye --nano-raw ${reads} \
      --asm-coverage 40 \
      --min-overlap 1000 \
      -g 100000 \
      -o assembly \
      --plasmids \
      -t ${task.cpus}
    """
    // TODO: add genome size est
}

process FAKE_ASSEMBLY {
    tag "${pair_id}" 

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("assembly"), emit: assembly

    script:
    """
    mkdir assembly
    gunzip -c ${reads} > assembly/${pair_id}_contigs.fna
    """
}

process PLASMID_ASSEMBLY{
    errorStrategy 'ignore'
    container "staphb/spades:latest"

    tag "${pair_id}"
    label 'big_mem'

    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'move', pattern: "plasmids/*"

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("plasmids"), emit: plasmids
    path("versions.yml"), emit: versions

    script:
    def single = reads instanceof Path

    def input = !single ? "-1 '${reads[0]}' -2 '${reads[1]}'" : "-s '${reads}'"
    def memory = "${task.memory}".replaceAll("[^0-9]","") 
    """
    spades.py -o tmp ${input} \
        -m ${memory} \
        -t ${task.cpus} \
        --plasmid
    mkdir plasmids
    mv tmp/contigs.fasta "plasmids/${pair_id}_contigs.fna"
    mv tmp/scaffolds.fasta "plasmids/${pair_id}_scaffolds.fna" || true
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spades: \$(spades.py -v | grep -oP '\\d{1,2}.\\d{1,2}.\\d{1,2}';)
    END_VERSIONS    
    """
}

workflow assembly_plasmids {
    take: reads
    main:
    // Plasmid assembly
    PLASMID_ASSEMBLY(reads)

    emit:
    versions = PLASMID_ASSEMBLY.out.versions
}

workflow assembly {    
    take: reads
    main:

    ch_versions = Channel.empty()

    if (params.assemblies) {
        reads
          .map{ tuple(it[0].replace(/\.fa/, "")
                           .replaceAll(/\./, "_"), it[1]) 
           }
          .set{ reads_clean_name }


        FAKE_ASSEMBLY(reads_clean_name)
        assemblies = FAKE_ASSEMBLY.out.assembly
    } else {
        // Shovil assembly
        ASSEMBLY(reads)
        ASSEMBLY.out.assembly.set{assemblies}
        ch_versions = ch_versions.mix(
            ASSEMBLY.out.versions.first()
        )
    }
    emit:
        assemblies = assemblies
        versions = ch_versions

}
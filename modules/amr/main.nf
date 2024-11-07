
process AMR_FINDER_UPDATE_DB {
    container null
    conda "bioconda::ncbi-amrfinderplus"

    cpus 1

    output:
    path("amrdb/*.1"), emit: amr_db

    script:
    """
    mkdir amrdb
    amrfinder_update --database amrdb
    """
}

process AMR_FINDER {
    container null
    conda "bioconda::ncbi-amrfinderplus"
    publishDir "${params.OUTDIR}/${params.RUNNAME}/${id}", mode: 'copy'

    cpus 4

    input:
    tuple val(id), path(assembly) 
    path(amr_db)

    output:
    path("*_amr_hits"), emit: amr
    path("versions.yml"), emit: versions
    script:
    def outf = assembly.baseName + "_amr_hits"
    """
    amrfinder \\
      --nucleotide $assembly/*.fna \\
      --database $amr_db \\
      --plus > $outf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinder: \$(amrfinder --version;)
    END_VERSIONS    
    """
}

process RESFINDER {

  container 'genomicepidemiology/resfinder:latest'
  publishdir "${params.OUTDIR}/${params.RUNNAME}/${id}", mode: 'copy'

  cpus 4

  input:
  tuple val(id), path(assembly), val(species)

  output:
  path("*_resfinder"), emit: amr
  path("versions.yml"), emit: versions

  script:
  def outf = assembly.baseName + "_resfinder"
  """
  resfinder -ifa $assembly -o $outf -s $species --acquired --point
  cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        refinder: \$(amrfinder --version;)
  END_VERSIONS
  """

}

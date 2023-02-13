// PARAMETERS
params.CONTAINER = "theoaphidian/wgs-illumina"

params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"
params.runName = "run01"

process READ_SAMPLESHEET {

    publishDir "${outputDir}", mode: 'copy'

    input:
    path(samplesheet) // input samplesheet
    path(outputDir) // where to save resulting tsv

    output:
    path "samplesheet.tsv", emit: samplesheetAbsolute
    
    script:
    """
    read_samplesheet.py -s "${params.samplesheet}" -i "${params.sampleName}" \
    -f "${params.fw_reads}" -r "${params.rv_reads}" -n "${params.runName}"
    """
}


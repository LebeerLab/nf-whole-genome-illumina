// PARAMETERS
params.CONTAINER = "theoaphidian/wgs-illumina"

params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"
params.runName = "run01"
params.single_end = false

process READ_SAMPLESHEET {

    publishDir "${outputDir}", mode: 'copy'

    input:
    path(samplesheet) // input samplesheet
    path(outputDir) // where to save resulting tsv

    output:
    path "samplesheet.tsv", emit: samplesheetAbsolute
    
    script:
    def single_end = params.single_end ? "-se" : ""
    """
    read_samplesheet.py read -s "${params.samplesheet}" -i "${params.sampleName}" \
    -f "${params.fw_reads}" -r "${params.rv_reads}" -n "${params.runName}" ${single_end}
    """
}


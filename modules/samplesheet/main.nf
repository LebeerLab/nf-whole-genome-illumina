// PARAMETERS
params.CONTAINER = "theoaphidian/wgs-illumina"

params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.sampleDatabase = "./sampledb"
params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"
params.runName = "run01"

process READ_SAMPLESHEET {

    publishDir "${outputDir}", mode: 'copy'

    input:
    path(samplesheet) // input samplesheet
    val(isCorrected) // correct samplesheet?
    path(outputDir) // where to save resulting tsv

    output:
    path "*.tsv", emit: samplesheetAbsolute
    
    script:
    def addToDB = isCorrected ? "${params.runName} --corrected" : "${params.runName}"
    """
    read_samplesheet.py "${params.samplesheet}" "${params.sampleName}" \
    "${params.fw_reads}" "${params.rv_reads}" "${params.sampleDatabase}" \
    "${addToDB}"
    """
}
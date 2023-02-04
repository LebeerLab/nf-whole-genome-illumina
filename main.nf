params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.outdir = "results"
params.sampleDatabase = "./sampledb"
params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"
params.runName = "run01"

params.debug = false

params.skip_fastp = false
params.truncLen = 0
params.trimLeft = 0
params.trimRight = 0
params.minLen = 50
params.maxN = 2

include { FASTP; MULTIQC } from './modules/qc' addParams(OUTPUT: "${params.outdir}")

def helpMessage() {
    log.info"""
     Name: nf-whole-genome-illumina
     Author: LAMB (UAntwerp)
    =========================================
    Required arguments:
      --samplesheet             Path to the samplesheet containing the reads to be assembled. Default: ${params.samplesheet}
      --sampleName              Column with the name of the sample. Default: ${params.sampleName}
      --fw_reads                Column name of the forward read. Default: ${params.fw_reads}
      --rv_reads               Column name of the reverse read. Default: ${params.rv_reads}
      --sampleDatabase          Path to a tsv file to store the location of the result and input. Default: ${params.sampleDatabase}
      --runName                 Name of the run. Default: ${params.runName}

    Optional arguments:

      --help  --h               Shows this help page
      --debug                   Run on a small subset of samples, for debugging purposes.
      --outdir                  The output directory where the results will be saved. Defaults to ${params.outdir}
      
      --truncLen                Truncation length used by fastp. Default = ${params.truncLen}
      --trimLeft --trimRight    Trimming on left or right side of reads by fastp. Default = ${params.trimLeft}
      --minLen                  Minimum length of reads kept by fastp. Default = ${params.minLen}
      --maxN                    Maximum amount of uncalled bases N to be kept by fastp. Default = ${params.maxN}

    Usage example:
        nextflow run main.nf --samplesheet '/path/to/samplesheet'
    """.stripIndent()
}

def paramsUsed() {
    log.info"""
    N F - W G S - I L L U M I N A
    =========================================
    samplesheet:      ${params.samplesheet}
    sampleName:       ${params.sampleName}
    runName:          ${params.runName}
    fw_reads:         ${params.fw_reads}
    rv_reads:         ${params.rv_reads}
    outdir:           ${params.outdir}
    sampleDatabase:   ${params.sampleDatabase}

    """.stripIndent()
}

if (params.help  || params.h){
    helpMessage()
    exit 0
}

process READ_SAMPLESHEET {

    publishDir "${params.sampleDatabase}", mode: 'copy'

    input:
    path(samplesheet) 

    output:
    path("samplesheet.tsv") 

    script:
    """
    read_samplesheet.py "${params.samplesheet}" "${params.sampleName}" \
    "${params.fw_reads}" "${params.rv_reads}" "${params.sampleDatabase}" "${params.runName}"
    """

}

workflow {
    
    paramsUsed()
    
    // Read samplesheet: find and update paths to reads (externalize from nf?)
    READ_SAMPLESHEET(params.samplesheet)
        .set{ samplesheetAbsolute }

    // Extract reads from samplesheet
    samplesheetAbsolute
        .splitCsv(header: true, sep: "\t") 
        .map {row -> tuple(row.ID, tuple(file(row.fw_reads), file(row.rv_reads)))}
        .set{reads_ch}

    if (!params.skip_fastp){

        //Filter and trim using fastp
        FASTP(reads_ch)

        FASTP.out.filteredReads
            .ifEmpty { error "No reads to filter"}
            .set { filteredReads }

        FASTP.out.fastp
            .collect()
            .set{fastp}

        MULTIQC(fastp)
    } else {
        filteredReads = reads_ch
    }

}

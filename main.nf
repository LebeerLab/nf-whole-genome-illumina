params.samplesheet = "${projectDir}/data/samplesheet.tsv"
params.outdir = "results"
params.sampleDatabase = "sampledb.tsv"
params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rev_reads = "rv_reads"
params.runName = "run01"

params.debug = false

params.skip_fastp = false
params.truncLen = 0
params.trimLeft = 0
params.trimRight = 0
params.minLen = 50
params.maxN = 2

def helpMessage() {
    log.info"""
     Name: nf-whole-genome-illumina
     Author: LAMB (UAntwerp)
    =========================================
    Required arguments:
      --samplesheet             Path to the samplesheet containing the reads to be assembled. Default: ${params.samplesheet}
      --sampleName              Column with the name of the sample. Default: ${params.sampleName}
      --fw_reads                Column name of the forward read. Default: ${params.fw_reads}
      --rev_reads               Column name of the reverse read. Default: ${params.rev_reads}
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
    "${params.fw_reads}" "${params.rev_reads}" "${params.sampleDatabase}" "${params.runName}"
    """

}

workflow {

    paramsUsed()

    READ_SAMPLESHEET(params.samplesheet)

    // Read samplesheet: find and update paths to reads (externalize from nf?)


    // Collect all fastq files
    // Channel
    //     .fromFilePairs(params.reads, size: params.pairedEnd ? 2 : 1)
    //     .ifEmpty { error "Could not find any reads matching the pattern ${params.reads}"}
    //     .take( params.debug ? 3 : -1 )
    //     //remove 'empty' samples
    //     .branch {
    //         success : params.pairedEnd ? it[1][1].countFastq() >= params.min_reads &&  it[1][0].countFastq() >= params.min_reads : it[1][0].countFastq() >= params.min_reads 
    //         failed : params.pairedEnd ? it[1][1].countFastq() < params.min_reads &&  it[1][0].countFastq() < params.min_reads : it[1][0].countFastq() < params.min_reads
    //     }
    //     .set { reads }

    // reads.failed.subscribe { println "Sample ${it[0]} did not meet minimum reads requirement of ${params.min_reads} reads and is excluded."}

    // if (!params.skip_fastp){

    //     // Filter and trim using fastp
    //     FASTP(reads.success)

    //     FASTP.out.filteredReads
    //         .ifEmpty { error "No reads to filter"}
    //         .set { filteredReads }

    //     FASTP.out.fastp
    //         .collect()
    //         .set{fastp}

    //     MULTIQC(fastp)
    // } else {
    //     filteredReads = reads.success
    // }

}

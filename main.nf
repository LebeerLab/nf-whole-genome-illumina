params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.outdir = "results"

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

//===== INCLUDE MODULES ==========================================
include { FASTP; MULTIQC } from './modules/qc' addParams(OUTPUT: "${params.outdir}")
include { READ_SAMPLESHEET } from './modules/samplesheet' addParams(
    sampleName : "${params.sampleName}", fw_reads : "${params.fw_reads}", 
    rv_reads : "${params.rv_reads}", runName : "${params.runName}")

// ================================================================
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
    """.stripIndent()
}

if (params.help  || params.h){
    helpMessage()
    exit 0
}

process ASSEMBLY {
    container "staphb/shovill:latest"

    tag "${pair_id}" 

    publishDir "${params.outdir}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("assembly")

    script:
    def single = reads instanceof Path

    def input = !single ? "--R1 '${reads[0]}' --R2 '${reads[1]}'" : "--R1 '${reads}'"
    """
    shovill --outdir assembly ${input} --ram ${task.memory} --cpus ${task.cpus}
    """
}

process CHECKM {
    container "nanozoo/checkm:latest"

    tag "${pair_id}" 

    publishDir "${params.outdir}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    tuple val(pair_id), path("results.tsv")

    script:
    """ 
    checkm lineage_wf -t ${task.cpus} ${assembly} lin --reduced_tree
    checkm qa lin/lineage.ms lin -t ${task.cpus} -o 2 --tab_table -f results.tsv
    """
}

workflow {
    
    paramsUsed()
    
    // Read samplesheet: find and update paths to reads (externalize from nf?)

    READ_SAMPLESHEET(params.samplesheet, file(params.samplesheet).getParent()) 

    // Extract reads from samplesheet
    READ_SAMPLESHEET.out.samplesheetAbsolute
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

    // Shovil assembly
    ASSEMBLY(filteredReads)
        .set{ assembly_ch }
    
    CHECKM(assembly_ch)

}

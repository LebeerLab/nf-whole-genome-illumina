params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.outdir = "results"

params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"
params.runName = "run01"
params.gtdb_db = null

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
      --gtdb_db                 Path to gtdb reference database. Default = ${params.gtdb_db}
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
    gtdb_db:          ${params.gtdb_db}
    """.stripIndent()
}

if (params.help  || params.h){
    helpMessage()
    exit 0
}

process ASSEMBLY {
    container "staphb/shovill:latest"

    tag "${pair_id}" 

    publishDir "${params.outdir}/${pair_id}/${params.runName}", mode: 'copy'

    input:
    tuple val(pair_id), path(reads)

    output:
    tuple val(pair_id), path("assembly")

    script:
    def single = reads instanceof Path

    def input = !single ? "--R1 '${reads[0]}' --R2 '${reads[1]}'" : "--R1 '${reads}'"
    """
    shovill --outdir assembly ${input} --ram ${task.memory} --cpus ${task.cpus}
    mv assembly/contigs.fa "assembly/${pair_id}_contigs.fna"
    """
}

process CHECKM {
    container "nanozoo/checkm:latest"

    //errorStrategy "ignore"
    tag "${pair_id}" 

    publishDir "${params.outdir}/${pair_id}/${params.runName}", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    tuple val(pair_id), path("checkm_results.tsv")

    script:
    """ 
    checkm lineage_wf ${assembly} lin --reduced_tree -t ${task.cpus}
    checkm qa lin/lineage.ms lin -t ${task.cpus} -o 2 --tab_table -f checkm_results.tsv
    """
}


process ANNOTATION {
    container "staphb/prokka"

    tag "${pair_id}"
    publishDir "${params.outdir}/${pair_id}/${params.runName}", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    tuple val(pair_id), path("annotation")
    script:
    """
    prokka "${assembly}/${pair_id}_contigs.fna" --outdir annotation --prefix "${pair_id}" \
    --locustag "${pair_id}" --compliant --cpus ${task.cpus} --quiet > prokka.out
    gzip -r "annotation"
    """
}

process DETECT_CHIMERS_CONTAMINATION {
    container "metashot/gunc"

    tag "${pair_id}"
    publishDir "${params.outdir}/${pair_id}/${params.runName}", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    tuple val(pair_id), path('gunc.out')

    script:
    when: params.gtdb_db != null
    """
    gunc -i assembly -r "${params.gtdb_db}"
    """
}

process MERGE_QC {
    container "metashot/gunc"

    tag "${pair_id}"
    publishDir "${params.outdir}/${pair_id}/${params.runName}", mode: 'copy'

    input:
    tuple val(pair_id), path(checkm_f), path(gunc_f)

    output:
    tuple val(pair_id), path(qc_f)

    script:
    when: params.gtdb_db != null
    """
    gunc merge_checkm --gunc_file gunc_f --checkm_file checkm_f --out_dir qc_f
    """
}

process CLASSIFICATION {
    container "theoaphidian/gtdbtk-entry"
    containerOptions "-v ${params.gtdb_db}:/refdata"

    tag "${pair_id}"
    publishDir "${params.outdir}/${pair_id}/${params.runName}", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)
    path(gtdb_db)

    output:
    tuple val(pair_id), path(output)
    script:
    when: params.gtdb_db != null
    """
    gtdbtk classify_wf --genome_dir assembly --prefix "${pair_id}_gtdbtk" --out_dir "output" --cpus $task.cpus --mash_db "mdb"
    """
}    


workflow {
    
    paramsUsed()
    
    // Read samplesheet: find and update paths to reads (externalize from nf?)
    if (!params.skip_samplesheet) {
    	READ_SAMPLESHEET(params.samplesheet, file(params.samplesheet).getParent()) 

    	// Extract reads from samplesheet
    	READ_SAMPLESHEET.out.samplesheetAbsolute
        	.splitCsv(header: true, sep: "\t") 
        	.map {row -> tuple(row.ID, tuple(file(row.fw_reads), file(row.rv_reads)))}
        	.set{reads_ch}
    } else {
	reads_ch = Channel.fromFilePairs(params.reads)
    }

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
    
    // QC  
    CHECKM(assembly_ch)
        .set{ checkm_ch }
    //DETECT_CHIMERS_CONTAMINATION(assembly_ch)
    //    .set { gunc_ch }
    //MERGE_QC(checkm_ch.join(gunc_ch))

    // Annotation genes
    ANNOTATION(assembly_ch)
    if (params.gtdb_db != null) {
    // Classification
    CLASSIFICATION(assembly_ch, file(params.gtdb_db))
    }
}

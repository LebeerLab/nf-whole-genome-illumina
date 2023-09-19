params.samplesheet = "${projectDir}/data/samplesheet.csv"
params.outdir = "results"


params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"
params.runName = "run01"

params.gtdb_db = null
params.mash_db = "mdb"
params.gunc_db = null
params.bakta_db = null

params.truncLen = 0
params.trimLeft = 0
params.trimRight = 0
params.minLen = 50
params.maxN = 2

params.depth = 150
params.minLenContig = 0
params.minCov = 2

params.skip_samplesheet = false
params.skip_fastp = false
params.skip_fastani = false
params.single_end = false
params.plasmids_only = false

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
      --rv_reads                Column name of the reverse read. Default: ${params.rv_reads}
      --runName                 Name of the run. Default: ${params.runName}
      --gtdb_db                 Path to gtdb reference database. Default = ${params.gtdb_db}
      --gunc_db                 Path to a gunc reference database. Default = ${params.gunc_db}
      --mash_db                 Path to the mash db made from the gtdb_db. Leave as is if you want the pipeline to generate this file (takes 45 minutes). Default = ${params.mash_db}
      --bakta_db                Path to the annotation db needed for bakta. Default = ${params.bakta_db}

    Optional arguments:

      --help  --h               Shows this help page
      --debug                   Run on a small subset of samples, for debugging purposes.
      --outdir                  The output directory where the results will be saved. Defaults to ${params.outdir}
      
      --single_end              If the input data is single end set this to true. Default = ${params.single_end}
      --truncLen                Truncation length used by fastp. Default = ${params.truncLen}
      --trimLeft --trimRight    Trimming on left or right side of reads by fastp. Default = ${params.trimLeft}
      --minLen                  Minimum length of reads kept by fastp. Default = ${params.minLen}
      --maxN                    Maximum amount of uncalled bases N to be kept by fastp. Default = ${params.maxN}
    
      --depth                   Subsample reads to this depth. Disable with --depth 0. Default = ${params.depth}
      --minLenContig            Minimum contig length. Set to 0 for automatic determination. Default = ${params.minLenContig}
      --minCov                  Minimumcontig coverage. Set to 0 for automatic determination. Default = ${params.minCov}
      --plasmids_only           Set to true if only the assemblies of the plasmids should be generated. Default = ${params.plasmids_only}
    
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
    gunc_db:          ${params.gunc_db}
    bakta_db:         ${params.bakta_db}
    """.stripIndent()
}

if (params.help  || params.h){
    helpMessage()
    exit 0
}

process ASSEMBLY {
    errorStrategy 'ignore'
    container "staphb/shovill:latest"

    tag "${pair_id}" 
    label 'big_mem'

    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'copy', pattern: "assembly/*"

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

process CHECKM {
    container "nanozoo/checkm:latest"

    tag "${pair_id}" 
    label 'big_mem'

    //publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(assembly)

    output:
    tuple val(pair_id), path("checkm_results.tsv"), emit: results
    path("versions.yml"), emit: versions

    script:
    """ 
    checkm lineage_wf ${assembly} lin --reduced_tree -t ${task.cpus}
    checkm qa lin/lineage.ms lin -t ${task.cpus} -o 2 --tab_table -f checkm_results.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkm: \$(checkm | grep -oP '\\d{1,2}.\\d{1,2}.\\d{1,2}';)
    END_VERSIONS    
    """
}

process ANNOTATION {
    container null
    //needs container locally: "oschwengers/bakta"
    tag "${pair_id}"

    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'copy', pattern: "annotation/*"

    input:
    tuple val(pair_id), path(assembly)

    output:
    tuple val(pair_id), path("annotation"), emit: annotation
    path('versions.yml'), emit: versions
    script:
    """
    bakta-docker.sh --db ${params.bakta_db} --output annotation \
    --prefix "${pair_id}" \
    --compliant --threads ${task.cpus} \
    "${assembly}/${pair_id}_contigs.fna"
    gzip annotation/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bakta: \$(docker run oschwengers/bakta --version | sed 'bakta //';)
    END_VERSIONS    
    """
}

process DETECT_CHIMERS_CONTAMINATION {
    container "metashot/gunc:1.0.5-1"

    tag "${pair_id}"

    input:
    tuple val(pair_id), path(assembly)
    path(guncdb)

    output:
    tuple val(pair_id), path("*.tsv"), emit: gunc
    path("versions.yml"), emit: versions

    script:
    """
    gunc run --input_dir assembly -r ${guncdb} \
        --file_suffix .fna \
        --threads ${task.cpus}        

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunc: \$(gunc -v ;)
    END_VERSIONS    
    """
}

process MERGE_QC {
    container "metashot/gunc:1.0.5-1"

    tag "${pair_id}"
    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'copy'

    input:
    tuple val(pair_id), path(checkm_f), path(gunc_f)

    output:
    path("qc/*.tsv")

    script:
    """
    mkdir qc
    gunc merge_checkm --gunc_file ${gunc_f} --checkm_file ${checkm_f} --out_dir qc
    """
}

process CLASSIFICATION {
    container "theoaphidian/gtdbtk-entry"
    containerOptions "-v ${params.gtdb_db}:/refdata"
    //conda "$projectDir/conda-env/gtdbtk.yaml"
    //publishDir "${params.outdir}/${params.runName}", mode: 'copy'
    label "big_mem"

    input:
    path(contigs)
    path(gtdb_db)
    path(mash_db)

    output:
    tuple path("*bac120.summary.tsv"), path("**bac120.ani_summary.tsv"), emit: classif
    path("versions.yml"), emit: versions
    //path("mash_db"), emit: mash_db

    script:
    def fastani = params.skip_fastani ? "--skip_ani_screen" : "--mash_db $mash_db"
    """
    gtdbtk classify_wf \
    --genome_dir . \
    --out_dir . \
    --cpus ${task.cpus} \
    --pplacer_cpus 1 \
    --scratch_dir tmp \
    $fastani
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gtdbtk: \$(gtdbtk -v | grep -oP '\\d{1,2}.\\d{1,2}.\\d{1,2}';)
    END_VERSIONS    
    """
}


process MERGE_CLASSIFICATION {

    publishDir "${params.outdir}/${params.runName}", mode: 'move'

    input:
    tuple path(bac_summ), path(bac_ani_summ)

    output:
    path("classification_summary.tsv")

    script:
    """
    merge_class.py
    """

}

process ANTISMASH {
    container null
    tag "${pair_id}"

    publishDir "${params.outdir}/${params.runName}/${pair_id}", mode: 'move', pattern: "antismash/*"

    input:
    tuple val(pair_id), path(annotation) 

    output:
    tuple val(pair_id), path("antismash/*"), emit: results
    path("versions.yml"), emit: versions
    script:
    """
    gunzip -c $annotation/*.gbff.gz > antismash.gbff
    run_antismash antismash.gbff . -c ${task.cpus} --genefinding-tool none

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(docker run antismash/standalone . --help | grep -oP '\\d{1,2}.\\d{1,2}.\\d{1,2}';)
    END_VERSIONS    
    """

}

process AMR_FINDER {
    container null
    conda "bioconda::ncbi-amrfinderplus"
    //publishDir "${params.outdir}/${params.runName}", mode: 'copy'

    input:
    path(assembly) 

    output:
    path("*_amr_hits"), emit: amr
    path("versions.yml"), emit: versions
    script:
    def outf = assembly.baseName + "_amr_hits"
    """
    amrfinder --nucleotide $assembly > $outf
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinder: \$(amrfinder --version;)
    END_VERSIONS    
    """
}

workflow read_samplesheet {
   main:

    READ_SAMPLESHEET(file(params.samplesheet, checkIfExists: true), file(params.samplesheet).getParent()) 

    // Extract reads from samplesheet
    READ_SAMPLESHEET.out.samplesheetAbsolute
        .splitCsv(header: true, sep: "\t") 
        .take( params.debug ? 2 : -1 )
	.set{reads_intermediate}

    if (params.single_end) {
        reads_intermediate
            .map {row -> tuple(row.ID, file(row.fw_reads))}
            .set{reads_ch}
    } else {
        reads_intermediate
            .map {row -> tuple(row.ID, tuple(file(row.fw_reads), file(row.rv_reads)))}
            .set{reads_ch}        
    
    }

    emit: 
        reads_ch
}

workflow filter_reads {
    take: reads
    main:

    ch_versions = Channel.empty()
    //Filter and trim using fastp
    FASTP(reads)
    ch_versions = ch_versions.mix(
        FASTP.out.versions.first()
    )

    FASTP.out.filteredReads
        .ifEmpty { error "No reads to filter"}
        .set { filteredReads }

    FASTP.out.fastp
        .collect()
        .set{fastp}

    MULTIQC(fastp)
    ch_versions = ch_versions.mix(
        MULTIQC.out.versions.first()
    )

    emit:
        reads = filteredReads
        versions = ch_versions
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
    // Shovil assembly
    ASSEMBLY(reads)

    ch_versions = ch_versions.mix(
        ASSEMBLY.out.versions.first()
    )

    contigs_ch = ASSEMBLY.out.assembly
                     .collect{it[1] + "/${it[0]}_contigs.fna"}
    
    // QC  
    CHECKM(ASSEMBLY.out.assembly)
    checkm_ch = CHECKM.out.results
    ch_versions = ch_versions.mix(
        CHECKM.out.versions.first()
    )

    if (params.gunc_db != null) {
        DETECT_CHIMERS_CONTAMINATION(ASSEMBLY.out.assembly, file(params.gunc_db))
        gunc_ch = DETECT_CHIMERS_CONTAMINATION.out.gunc
        ch_versions = ch_versions.mix(
            DETECT_CHIMERS_CONTAMINATION.out.versions.first()
        )
        MERGE_QC(checkm_ch.join(gunc_ch))
            .collectFile(
                keepHeader: true, 
                name: "qc_checkm_gunc.tsv", 
                storeDir: "${params.outdir}/${params.runName}")
    }
    if (params.bakta_db != null) {
        // Annotation genes
        ANNOTATION(ASSEMBLY.out.assembly)
        ANNOTATION.out.annotation
            // grab AMB.gbk.gz file from output
            .map { it -> [it[0], 
                        file(it[1] + "/${it[0]}.gbff.gz")
                        ] 
                }
        predicted_genes_ch = ANNOTATION.out.annotation
        ch_versions = ch_versions.mix(
            ANNOTATION.out.versions.first()
        )

        ANTISMASH(predicted_genes_ch)
        ch_versions = ch_versions.mix(
            ANTISMASH.out.versions.first()
        )
    }
    emit:
        contigs = contigs_ch
        versions = ch_versions
}

workflow classification {
    take: contigs
    main:
    if (params.gtdb_db != null) {
      
        ch_versions = Channel.empty()
        CLASSIFICATION(contigs, file(params.gtdb_db), file(params.mash_db))
        ch_versions = ch_versions.mix(
            CLASSIFICATION.out.versions.first()
        )
        MERGE_CLASSIFICATION(CLASSIFICATION.out.classif)
    }
    emit:
    versions = ch_versions
}

workflow {
    paramsUsed()
    ch_versions = Channel.empty()
    read_samplesheet()

    def execute_fastp = params.skip_fastp ? false : true
    if (execute_fastp) {
        filter_reads(read_samplesheet.out)
        reads = filter_reads.out.reads
        ch_versions = ch_versions.mix(
            filter_reads.out.versions
        )
    } else {
        reads = read_samplesheet.out
    }
    
    if (params.plasmids_only){
        assembly_plasmids(reads)
        ch_versions = ch_versions.mix(
            assembly_plasmids.out.versions
        )
    } else {
        assembly(reads)
        ch_versions = ch_versions.mix(
            assembly.out.versions
        )
        assembly_plasmids(reads)
        ch_versions = ch_versions.mix(
            assembly_plasmids.out.versions
        )
        classification(assembly.out.contigs)
        ch_versions = ch_versions.mix(
            classification.out.versions
        )
        AMR_FINDER(assembly.out.contigs)
        ch_versions = ch_versions.mix(
            AMR_FINDER.out.versions.first()
        )
        AMR_FINDER.out.amr
         .collectFile(name: 'amrfinder_results.tsv', storeDir: params.outdir, keepHeader: true)
    }

    ch_versions.unique().collectFile(name: 'software_versions.yml', storeDir: params.outdir)
}

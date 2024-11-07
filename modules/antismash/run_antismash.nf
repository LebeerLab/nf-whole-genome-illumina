params.fnadir = "${launchDir}/fnas/*.fna"


include { ANTISMASH } from './main.nf' addParams(OUTDIR: "${params.outdir}", RUNNAME: "${params.runName}")
include { ANNOTATION } from './../../main.nf'

workflow {

    fnadir = channel.fromPath(params.fnadir)
    fnas = fnadir.map { tuple(file(it).getBaseName(), it ) } 
    ANNOTATION(fnas)
    ANTISMASH(ANNOTATION.out.annotation)
}

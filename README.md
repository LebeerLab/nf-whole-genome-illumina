# nf-whole-genome-illumina

## Samplesheet input

The samplesheet is updated and the relative location of the fastq files can be deduced starting from the folder above the samplesheet. 
Therefore it is important to supply the samplesheet in the following format: `dir/containing/samplesheet/samplesheet.tsv`.
The pipeline accepts tsv or csv samplesheets, as it will convert csv sheets to tsv.

The samplesheet should contain atleast the three columns with the following names in the header:
params.sampleName = "ID"
params.fw_reads = "fw_reads"
params.rv_reads = "rv_reads"


## Reference Databases

This pipeline needs three separate reference databases for the whole pipeline to run:

* gtdb_db: Needed for classification using GTDBTk. [Needs to be downloaded and untarred](https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz)
* mash_db: The mash database constructed from the earlier mentioned GTDB database is stored. If not supplied, the pipeline will generate this by itself, but this will take aproximately 45 minutes. 
* gunc_db: Needed for the gunc quality check. Can be downloaded by `gunc download_db`, or [retrieved manually here](https://swifter.embl.de/~fullam/gunc/). Usethe progenomes database. It is smaller, faster and should produce similar results. 
* bakta_db: Needed for annotation using the [bakta](https://github.com/oschwengers/bakta) tool. Can be downloaded [here](https://zenodo.org/record/7669534).

## Dependencies

- Bakta requires a local docker version via

`docker pull oschwengers/bakta`

due to the the way it is set up.
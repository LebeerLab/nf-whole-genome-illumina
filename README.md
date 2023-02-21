# nf-whole-genome-illumina

## Reference Databases

This pipeline needs three separate reference databases for the whole pipeline to run:

* gtdb_db: Needed for classification using GTDBTk. [Needs to be downloaded and untarred](https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_v2_data.tar.gz)
* mash_db: The mash database constructed from the earlier mentioned GTDB database is stored. If not supplied, the pipeline will generate this by itself, but this will take aproximately 45 minutes. 
* gunc_db: Needed for the gunc quality check. Can be downloaded by `gunc download_db`, or [retrieved manually here](https://swifter.embl.de/~fullam/gunc/). Usethe progenomes database. It is smaller, faster and should produce similar results. 

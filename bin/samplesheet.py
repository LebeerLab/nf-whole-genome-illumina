import os
from pathlib import Path
import re
import pandas as pd
from utils import generate_uqid

CORR_SAMPLESHEET = "samplesheet.tsv"
DEF_SAMPLE_ID = "ID"
DEF_FW_READS = "fw_reads"
DEF_RV_READS = "rv_reads"
DEF_ASSEMBLY = "assembly"
DEF_RUN = "run01"

def error_message(message):
    red = '\033[91m'
    end = '\033[0m'
    return f"{red}{message}{end}"


class SampleSheet:
    def __init__(
        self,
        samplesheet=CORR_SAMPLESHEET,
        sample_col=DEF_SAMPLE_ID,
        fw_col=DEF_FW_READS,
        rev_col=DEF_RV_READS,
        sample_db_dir=None,
        run_id=None,
        paired_end=True
    ) -> None:
        # Validation input
        self.paired_end = paired_end
        if sample_col == None:
            sample_col = DEF_SAMPLE_ID
        
        ## Assume if fw_col is given, that this is on purpose (eg single reads)
        if fw_col == None:
            fw_col = DEF_FW_READS
        if paired_end and rev_col == None:
            rev_col = DEF_RV_READS
        
        if run_id == None:
            run_id = DEF_RUN

        self.sample_col = sample_col
        self.fw_col = fw_col
        self.rev_col = rev_col
        self.run_id = run_id
        self.root_dir = os.getcwd()
        self.paired_end = paired_end

        if sample_db_dir is not None:
            self.corrected_sheet = True
            self.sample_db_samplesheet = os.path.join(sample_db_dir, "sampledb.tsv")
            self.filename = os.path.join(os.path.dirname(samplesheet), CORR_SAMPLESHEET)
            samples_db_dir_path = Path(sample_db_dir)

            if not samples_db_dir_path.exists():
                os.makedirs(samples_db_dir_path, exist_ok=True)
        else:
            self.filename = samplesheet
            self.corrected_sheet = False

    def read_samplesheet(self):

        try:
            if self.filename.endswith(".tsv"):
                smpsh = pd.read_table(self.filename)
                
            ## Asumption: it is a csv, could implement xlsx with pd.read_xlsx,
            ## But then need to supply pyopenxlsx dependency in docker...
            else:
                smpsh = pd.read_csv(self.filename)
            ## Assertion of sample cols
            assert self.fw_col in smpsh.columns, f"Forward column {self.fw_column} not found in file with header {smpsh.columns.values}."
            if self.paired_end:
                assert self.rev_col in smpsh.columns, f"Reverse column {self.rev_col} not found in file with header {smpsh.columns.values}."
            self.content = smpsh
        except AssertionError as e:
            print(error_message(f"Error during reading samplesheet:\n{e}"))
            exit(1)

    def _build_read_paths(self, absolute=False):
        # Update paths to rel paths.
        if self.paired_end:
            self.content[[self.fw_col, self.rev_col]] = self.content[
                [self.fw_col, self.rev_col]
            ].apply(lambda x: self._fetch_filepath(x, absolute))
        else:
            self.content[[self.fw_col]] = self.content[
                    [self.fw_col]
                ].apply(lambda x: self._fetch_filepath(x, absolute))
                

    def _build_assembly_paths(self):
        # Update paths to absolute paths.
        self.content[[DEF_ASSEMBLY]] = self.content[[DEF_ASSEMBLY]].apply(
            lambda x: self._fetch_filepath(x, absolute=True, is_assembly=True)
        )

    def update_samplesheet(self):

        self._build_read_paths()
        # Add run_name as column
        self.content["run_id"] = self.run_id
        # Generate uqid
        uqid_columns = ["run_id", self.fw_col]
        if self.paired_end:
            uqid_columns.append(self.rev_col)
        self.content["uqid"] = (
            self.content[uqid_columns]
            .astype(str)
            .agg("".join, axis=1)
            .map(generate_uqid)
        )

        # Rename ID, rv_read, fw_read cols to fixed names
        colnames = {
                self.sample_col: DEF_SAMPLE_ID,
                self.fw_col: DEF_FW_READS,
            }
        if self.paired_end:
            colnames[self.rev_col] = DEF_RV_READS

        self.content.rename(
            columns=colnames
        )

        # Add resulting assembly path as column
        # RESULTS / RUN / ID/ ASSEMBLY / ID_contigs.fna
        root_contig = self.root_dir + "/results/"
        self.content[DEF_ASSEMBLY] = [
            f"{root_contig}{run}/{id}/assembly/{id}_contigs.fna"
            for id, run in zip(self.content[DEF_SAMPLE_ID], self.content["run_id"])
        ]

    def write_samplesheet(self):
        self.content.to_csv(CORR_SAMPLESHEET, sep="\t", index=False)

    def merge_summaries(self):
        
        # Assumption: data will originate from a single run most of the time
        # Therefore: load summary data in dict of pd objects and fetch results from those

        folder_col = "run_folder"
        checkm_col = "checkm_gunc"
        gtdb_col = "gtdb_bac120"
        
        col_by_name = {}
        col_by_name[checkm_col] = "/qc_checkm_gunc.tsv"
        col_by_name[gtdb_col] = "/gtdbtk.bac120.ani_summary.tsv"
        

        summary_data = {}
        summary_data[checkm_col] = {}
        summary_data[gtdb_col] = {}
        ## filename : pd.DataFrame

        new_data = self.content.copy() 
        # ../PROJ/RESULTS/RUN/SAMPLE/ASSEMBLY/contig.fna
        new_data[folder_col] = new_data[DEF_ASSEMBLY]\
            .str.split("/").str[:-3]\
                .str.join("/")
        
        df_run_folders = new_data.groupby(folder_col).first()
        for colname in [checkm_col, gtdb_col]:
            df_run_folders[colname] = df_run_folders.index.astype(str) + col_by_name[colname]
        
        for idx, row in df_run_folders.iterrows():
            summary_data[checkm_col][row[checkm_col]] = pd.read_table(row[checkm_col])
            summary_data[gtdb_col][row[gtdb_col]] = pd.read_table(row[gtdb_col])

        # Join back to copy of original data
        cols_to_keep = df_run_folders.columns.difference(new_data.columns)
        new_data.index = new_data[folder_col]
        new_data = pd.merge(new_data, df_run_folders[cols_to_keep], left_index=True, right_index=True, how="left")
        
        all_data = {}
        all_data[checkm_col] = []
        all_data[gtdb_col] = []

        def _fetch_result_from_df(row, colname_df):
            sample_id = row[DEF_SAMPLE_ID] + "_contigs"
            if colname_df == gtdb_col:
                df = summary_data[colname_df][row[gtdb_col]]
                return df.loc[df["user_genome"] == sample_id,]
            elif colname_df == checkm_col:
                df = summary_data[colname_df][row[checkm_col]]
                return df.loc[df["genome"] == sample_id,]
        
        
        for idx, row in new_data.iterrows():
            all_data[checkm_col].append(_fetch_result_from_df(row, checkm_col))
            all_data[gtdb_col].append(_fetch_result_from_df(row, gtdb_col))
        
        df_checkm = pd.concat(all_data[checkm_col])
        df_gtdb = pd.concat(all_data[gtdb_col])
        self.content = pd.concat([self.content, df_checkm, df_gtdb], axis=1)
        

    def update_sampledb(self):
        # Add absolute root to reads and assembly
        self._build_read_paths(absolute=True)
        self._build_assembly_paths()
        self.merge_summaries()
        output = self.content

        if Path(self.sample_db_samplesheet).exists():
            db_content = pd.read_table(self.sample_db_samplesheet)
            output = pd.concat(
                [
                    db_content,
                    self.content[~self.content["uqid"].isin(db_content["uqid"])],
                ]
            )

        output.to_csv(self.sample_db_samplesheet, sep="\t", index=False)

    def _fetch_filepath(self, inp_files, absolute=False, is_assembly=False):
        def simplify_samplenames(filenames: pd.Series) -> pd.Series:
            return filenames.apply(
                lambda x: x.lower().split(".")[0].split("/")[-1].strip()
            )

        def simplify_samplename(filename: str) -> str:
            return filename.lower().split(".")[0].strip()

        def find_matching_file(filename, rootdir, absolute=False):
            for root, _, files in os.walk(rootdir):
                for file in files:
                    if simplify_samplename(file) == filename:
                        if absolute:
                            root = os.path.join(self.root_dir, root)
                        return os.path.join(root, file)
            raise FileNotFoundError(
                f"Could not locate the file specified as {filename}"
            )

        if is_assembly:
            rt = os.path.join(self.root_dir, "results")
        else:
            rt = os.path.dirname(self.filename)
        simple_names = simplify_samplenames(inp_files)

        return simple_names.apply(lambda x: find_matching_file(x, rt, absolute))

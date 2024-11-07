#!/bin/python

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
DEF_ASS_METH = "assembly_method"
DEF_RUN = "run01"


def fatal_error_message(message):
    red = "\033[91m"
    end = "\033[0m"
    print(f"{red}{message}{end}")
    exit(1)

def clean_ambi(ambi:str):

    ambi = str(ambi).upper()
    try:
        amb_type, amb_no, amb_ext = re.findall("(^AMB[-]?[A-Z])[-]?(\d+)(.*)$", ambi)[0]
    except IndexError as e:
        #raise NameError(f"Could not identify the AMB code {ambi}")
        return ambi
    amb_type = amb_type.replace("-", "")
    amb_no = amb_no.zfill(4)
    if len(amb_ext) >0:
        amb_ext = f"-{amb_ext.replace('-','')}"
    return f"{amb_type}-{amb_no}{amb_ext}"


class SampleSheet:
    def __init__(
        self,
        samplesheet=CORR_SAMPLESHEET,
        sample_col=DEF_SAMPLE_ID,
        fw_col=DEF_FW_READS,
        rev_col=DEF_RV_READS,
        sample_db_dir=None,
        run_id=None,
        paired_end=True,
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

        # Write mode

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
            
            assert (
                self.fw_col in smpsh.columns
            ), f"Forward column {self.fw_col} not found in file with header {smpsh.columns.values}."
            if self.paired_end:
                assert (
                    self.rev_col in smpsh.columns
                ), f"Reverse column {self.rev_col} not found in file with header {smpsh.columns.values}."
            self.content = smpsh
        except AssertionError as e:
            fatal_error_message(f"Error during reading samplesheet:\n{e}")

    def _build_read_paths(self, absolute=False):
        # Update paths to rel paths.
        if self.paired_end:
            self.content[[self.fw_col, self.rev_col]] = self.content[
                [self.fw_col, self.rev_col]
            ].apply(lambda x: self._fetch_filepath(x, absolute))
        else:
            self.content[[self.fw_col]] = self.content[[self.fw_col]].apply(
                lambda x: self._fetch_filepath(x, absolute)
            )

    def _build_assembly_paths(self):
        # Update paths to absolute paths.
        self.content[[DEF_ASSEMBLY]] = self.content[[DEF_ASSEMBLY]].apply(
            lambda x: self._fetch_filepath(x, absolute=True, is_assembly=True)
        )
        
    def update_samplesheet(self):

        self._build_read_paths()
        print(self.content)
        # Add run_name as column
        self.content["run_id"] = self.run_id

        # Rename ID, rv_read, fw_read cols to fixed names
        colnames = {
            self.sample_col: DEF_SAMPLE_ID,
            self.fw_col: DEF_FW_READS,
        }
        if self.paired_end:
            colnames[self.rev_col] = DEF_RV_READS

        self.content.rename(columns=colnames)
        
        # Clean AMB identifiers where possible
        self.content[DEF_SAMPLE_ID] = self.content[DEF_SAMPLE_ID].apply(lambda x: clean_ambi(x))

        # Add resulting assembly path as column
        # RESULTS / RUN / ID/ ASSEMBLY / ID_contigs.fna
        root_contig = self.root_dir + "/results/"
        self.content[DEF_ASSEMBLY] = [
            f"{root_contig}{run}/{id}/assembly/{id}_contigs.fna.gz"
            for id, run in zip(self.content[DEF_SAMPLE_ID], self.content["run_id"])
        ]

    def write_samplesheet(self):
        self.content.to_csv(CORR_SAMPLESHEET, sep="\t", index=False)


    def merge_summaries(self):

        def set_sample_index(data):
            data.index = pd.MultiIndex.from_arrays(
                data[[DEF_SAMPLE_ID, "run_id"]].values.T,
                names=[DEF_SAMPLE_ID, "run_id"]
            )
            return data

        # Assumption: data will originate from a single run, with a single summary for the qc
        # and a single summary for the classification
        # Therefore: load both summary data in dict of pd objects and fetch results from those

        checkm_col = "checkm_gunc"
        gtdb_col = "gtdb_bac120"

        summary_data = {}
        summary_data[checkm_col] = {}
        summary_data[gtdb_col] = {}
        ## filename : pd.DataFrame

        new_data = self.content.copy()
        # ../PROJ/RESULTS/RUN/SAMPLE/ASSEMBLY/contig.fna
        assembly_roots = new_data[DEF_ASSEMBLY].str.split("/").str[:-3].str.join("/").unique()
        assert len(assembly_roots) == 1, f"More than one run supplied in same samplesheet!"

        checkm_file = f"{assembly_roots[0]}/qc_checkm_gunc.tsv"
        # Depends on the output obtained from classification...
        for filename in os.listdir(self.root_dir + "/results/" + self.run_id):
            if filename.endswith("summary.tsv"):
                classification_tsv = filename
        try:
            class_file = f"{assembly_roots[0]}/{classification_tsv}"
        except:
            print("No classification summary found for {self.filename}, aborting addition...")
            exit()
        checkm_data = pd.read_table(checkm_file)
        checkm_data[DEF_SAMPLE_ID] = checkm_data["genome"].str.replace("_contigs", "")
        checkm_data = checkm_data.drop(columns="genome")
        
        class_data = pd.read_table(class_file)
        class_data[DEF_SAMPLE_ID] = class_data["user_genome"].str.replace("_contigs", "")
        class_data = class_data.drop(columns="user_genome")

        new_data = set_sample_index(new_data)
        checkm_data["run_id"] = self.run_id
        checkm_data = set_sample_index(checkm_data)
        checkm_data = checkm_data.drop(columns=["run_id", DEF_SAMPLE_ID])
        class_data["run_id"] = self.run_id
        class_data = set_sample_index(class_data)
        class_data = class_data.drop(columns=["run_id", DEF_SAMPLE_ID])
        self.content = pd.concat([new_data, checkm_data, class_data], axis=1)

    def update_sampledb(self):
        # Add absolute root to reads and assembly
        self._build_read_paths(absolute=True)
        self._build_assembly_paths()

        # Generate uqid
        uqid_columns = ["run_id", DEF_FW_READS]

        self.content["uqid"] = (
            self.content[uqid_columns]
            .astype(str)
            .agg("".join, axis=1)
            .map(generate_uqid)
        )

        def extract_assembly_method():
            wd = os.getcwd()
            sr_assembly = "shovill 1.1.0"
            lr_assembly = "Flye 2.8.1-b1676"
            if re.search("nanopore",wd,re.IGNORECASE):
                return lr_assembly
            return sr_assembly

        def extract_child_folder_name(read_paths: pd.Series, parent_folder):
            try:
                return (
                    read_paths.str.split("/" + parent_folder + "/")
                    .str[-1]
                    .str.split("/")
                    .str[0]
                )
            except Exception as e:
                fatal_error_message(f"Parent folder '{parent_folder}' not found. {e}")

        # Extract platform
        self.content["platform"] = extract_child_folder_name(
            self.content[DEF_FW_READS], "seqdata"
        )

        # Extract project name
        self.content["project"] = extract_child_folder_name(
            self.content[DEF_FW_READS], "projects"
        )

        # Extract assembly method
        self.content[[DEF_ASS_METH]] = extract_assembly_method()

        # Add QC and class info
        try:
            self.merge_summaries()
        except FileNotFoundError as e:
            fatal_error_message(f"Results folder not found! {e}")

        # Remove following fields
        cols_to_remove = [
            "genome",
            "GUNC.n_contigs",
            "GUNC.n_genes_called",
            "GUNC.n_genes_mapped",
            "GUNC.divergence_level",
            "pass.MIMAG_medium",
            "pass.MIMAG_high",
            "user_genome",
        ]
        self.content = self.content.drop(columns=cols_to_remove, errors="ignore")

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
                lambda x: str(x).lower().split(".")[0].split("/")[-1].strip()
            )

        def simplify_samplename(filename: str) -> str:
            return filename.lower().split(".")[0].strip()

        def find_matching_file(filename, rootdir, absolute=False):
            for root, _, files in os.walk(rootdir):
                for file in files:
                    if simplify_samplename(file) == filename:
                        if absolute:
                            root = os.path.join(self.root_dir, root)
                        if is_assembly and file.endswith('.fna.gz'):
                           continue
                        return os.path.join(root, file)

            raise FileNotFoundError(
                f"Could not locate the file specified as {filename}; settings = {rootdir} absolute:{absolute}, is_assembly:{is_assembly}"
            )

        if is_assembly:
            rt = os.path.join(self.root_dir, "results")
        else:
            rt = os.path.dirname(self.filename)
        simple_names = simplify_samplenames(inp_files)

        return simple_names.apply(lambda x: find_matching_file(x, rt, absolute))

if __name__ == "__main__":
    import sys
    input = sys.argv[0]

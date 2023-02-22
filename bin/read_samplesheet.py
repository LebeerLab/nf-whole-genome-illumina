#!/usr/bin/env python3

import os
import argparse
from pathlib import Path
import pandas as pd

from utils import generate_uqid

CORR_SAMPLESHEET = "samplesheet.tsv"
DEF_SAMPLE_ID = "ID"
DEF_FW_READS = "fw_reads"
DEF_RV_READS = "rv_reads"
DEF_ASSEMBLY = "assembly"


class SampleSheet:
    def __init__(
        self,
        samplesheet=CORR_SAMPLESHEET,
        sample_col=DEF_SAMPLE_ID,
        fw_col=DEF_FW_READS,
        rev_col=DEF_RV_READS,
        sample_db_dir=None,
        run_id=None,
    ) -> None:
        self.sample_col = sample_col
        self.fw_col = fw_col
        self.rev_col = rev_col
        self.run_id = run_id
        self.root_dir = os.getcwd()

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

        if self.filename.endswith(".tsv"):
            smpsh = pd.read_table(self.filename)
            assert self.rev_col in smpsh.columns
        ## Asumption: it is a csv, could implement xlsx with pd.read_xlsx,
        ## But then need to supply pyopenxlsx dependency in docker...
        else:
            smpsh = pd.read_csv(self.filename)
        self.content = smpsh

    def _build_read_paths(self, absolute=False):
        # Update paths to absolute paths.
        self.content[[self.fw_col, self.rev_col, DEF_ASSEMBLY]] = self.content[
            [self.fw_col, self.rev_col, DEF_ASSEMBLY]
        ].apply(lambda x: self._fetch_filepath(x, absolute))


    def update_samplesheet(self):

        self._build_read_paths()
        # Add run_name as column
        self.content["run_id"] = self.run_id
        # Generate uqid
        self.content["uqid"] = (
            self.content[["run_id", self.fw_col, self.rev_col]]
            .sum(axis=1)
            .map(generate_uqid)
        )

        # Rename ID, rv_read, fw_read cols to fixed names
        self.content.rename(
            columns={
                self.sample_col: DEF_SAMPLE_ID,
                self.fw_col: DEF_FW_READS,
                self.rev_col: DEF_RV_READS,
            }
        )

        # Add resulting assembly path as column
        # RESULTS / ID / RUN / ASSEMBLY / ID_contigs.fna
        root_contig = os.path.dirname(self.filename).split("/data")[0] + "/results/"
        self.content[DEF_ASSEMBLY] = [
            f"{root_contig}{id}/{run}/assembly/{id}_contigs.fna"
            for id, run in zip(self.content[DEF_SAMPLE_ID], self.content["run_id"])
        ]

    def write_samplesheet(self):
        self.content.to_csv(CORR_SAMPLESHEET, sep="\t", index=False)

    def update_sampledb(self):
        # Add absolute root
        self._build_read_paths(absolute=True)

        if Path(self.sample_db_samplesheet).exists():
            db_content = pd.read_table(self.sample_db_samplesheet)
            db_content = db_content.append(
                self.content[~self.content[DEF_SAMPLE_ID].isin(db_content[DEF_SAMPLE_ID])],
                ignore_index=True,
            )

            db_content.to_csv(self.sample_db_samplesheet, sep="\t", index=False)
        else:
            self.content.to_csv(self.sample_db_samplesheet, sep="\t", index=False)

    def _fetch_filepath(self, read_files, absolute=False):
        def simplify_samplenames(filenames: pd.Series) -> pd.Series:
            return filenames.apply(
                lambda x: x.lower().split(".")[0].split("/")[-1].strip()
            )

        def simplify_samplename(filename: str) -> str:
            return filename.lower().split(".")[0].strip() 

        def find_matching_file(filename, rootdir, absolute=False):
            for root, _, files in os.walk(rootdir):
                for file in files:
                    if simplify_samplename(file) == simplify_samplename(filename):
                        if absolute:
                            absolute_root = os.path.join(self.root_dir, root)
                            return os.path.join(absolute_root, file)    
                        return os.path.join(root, file)
            raise FileNotFoundError(
                f"Could not locate the file specified as {filename}"
            )

        smpsh_dir = os.path.dirname(self.filename)
        reads_smp = simplify_samplenames(read_files)
        return reads_smp.apply(lambda x: find_matching_file(x, smpsh_dir, absolute))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="WGS Samplesheet Reader",
        description="Read WGS Samplesheets. Find samples listed and save to a masterfile.",
    )

    parser.add_argument("-s", "--samplesheet", default=CORR_SAMPLESHEET)
    parser.add_argument("-i", "--sample_column", default=DEF_SAMPLE_ID)
    parser.add_argument("-f", "--forward_column", default=DEF_FW_READS)
    parser.add_argument("-r", "--reverse_column", default=DEF_RV_READS)
    parser.add_argument("-d", "--sample_db_dir", default=None)
    parser.add_argument("-n", "--run_name", default=None)

    args = parser.parse_args()
    smpsh = SampleSheet(
        args.samplesheet,
        args.sample_column,
        args.forward_column,
        args.reverse_column,
        args.sample_db_dir,
        args.run_name,
    )

    smpsh.read_samplesheet()

    # Run these at the start of the pipeline
    if not smpsh.corrected_sheet:
        smpsh.update_samplesheet()
        smpsh.write_samplesheet()
    else:
        # Run this when finished
        smpsh.update_sampledb()

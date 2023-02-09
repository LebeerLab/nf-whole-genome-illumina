#!/usr/bin/env python3

import os
import argparse
from pathlib import Path
import pandas as pd

from utils import generate_uqid

CORR_SAMPLESHEET = "samplesheet.tsv"

class SampleSheet:
    def __init__(
        self,
        samplesheet,
        sample_col,
        fw_col,
        rev_col,
        sample_db_dir,
        run_id,
        corrected_sheet=False,
    ) -> None:
        self.sample_col = sample_col
        self.fw_col = fw_col
        self.rev_col = rev_col
        self.sample_db_samplesheet = os.path.join(sample_db_dir, "sampledb.tsv")
        self.run_id = run_id
        self.corrected_sheet = corrected_sheet
        if corrected_sheet:
            self.filename = os.path.join(os.path.dirname(samplesheet), CORR_SAMPLESHEET)
        else:
            self.filename = samplesheet
        
        samples_db_dir_path = Path(sample_db_dir)

        if not samples_db_dir_path.exists():
            os.makedirs(samples_db_dir_path, exist_ok=True)

    def read_samplesheet(self):

        if self.filename.endswith(".tsv"):
            smpsh = pd.read_table(self.filename)
            assert self.rev_col in smpsh.columns
        ## Asumption: it is a csv, could implement xlsx with pd.read_xlsx,
        ## But then need to supply pyopenxlsx dependency in docker...
        else:
            smpsh = pd.read_csv(self.filename)
        self.content = smpsh

    def update_samplesheet(self):
        # Update paths to absolute paths.
        self.content[[self.fw_col, self.rev_col]] = self.content[
            [self.fw_col, self.rev_col]
        ].apply(lambda x: self._fetch_filepath(x))

        # Add run_name as column
        self.content["run_id"] = self.run_id
        # Generate uqid
        self.content["uqid"] = (
            self.content[["run_id", self.fw_col, self.rev_col]]
            .sum(axis=1)
            .map(generate_uqid)
        )

        # TODO Add resulting assembly path as column
        # self.content["assembly"] = ""

        # Rename ID, rv_read, fw_read cols to fixed names
        self.content.rename(
            columns={
                self.sample_col: "ID",
                self.fw_col: "fw_reads",
                self.rev_col: "rv_reads",
            }
        )

    def write_samplesheet(self):
        self.content.to_csv(CORR_SAMPLESHEET, sep="\t", index=False)

    def update_sampledb(self):

        if Path(self.sample_db_samplesheet).exists():
            db_content = pd.read_csv(self.sample_db_samplesheet)
            db_content = db_content.append(
                self.content[~self.content["ID"].isin(db_content["ID"])],
                ignore_index=True,
            )

            db_content.to_csv(self.sample_db_samplesheet, sep="\t", index=False)
        else:
            self.content.to_csv(self.sample_db_samplesheet, sep="\t", index=False)

    def _fetch_filepath(self, read_files):
        def simplify_samplenames(filenames: pd.Series) -> pd.Series:
            return filenames.apply(
                lambda x: x.lower().split(".")[0].split("/")[-1].strip()
            )

        def simplify_samplename(filename: str) -> str:
            return filename.lower().split(".")[0].strip()

        def find_matching_file(filename, rootdir):
            for root, _, files in os.walk(rootdir):
                for file in files:
                    if simplify_samplename(file) == filename:
                        return os.path.join(root, file)
            raise FileNotFoundError(
                f"Could not locate the read specified as {filename}"
            )

        smpsh_dir = os.path.dirname(self.filename)
        reads_smp = simplify_samplenames(read_files)
        return reads_smp.apply(lambda x: find_matching_file(x, smpsh_dir))


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog="WGS Samplesheet Reader",
        description="Read WGS Samplesheets. Find samples listed and save to a masterfile.",
    )

    KNOWN_ARGS = ("samplesheet", "sample_column", "forward_column",
                "reverse_column", "sample_db_dir", "run_name", "corrected")
    parser.add_argument("samplesheet")
    parser.add_argument("sample_column")
    parser.add_argument("forward_column")
    parser.add_argument("reverse_column")
    parser.add_argument("sample_db_dir")
    parser.add_argument("run_name")
    parser.add_argument("-c", "--corrected", action="store_true")

    args = parser.parse_args()
    smpsh = SampleSheet(
        args.samplesheet,
        args.sample_column,
        args.forward_column,
        args.reverse_column,
        args.sample_db_dir,
        args.run_name,
        args.corrected,
    )

    smpsh.read_samplesheet()

    # Run these at the start of the pipeline
    if not smpsh.corrected_sheet:
        smpsh.update_samplesheet()
        smpsh.write_samplesheet()
    else:
        # Run this when finished
        smpsh.update_sampledb()

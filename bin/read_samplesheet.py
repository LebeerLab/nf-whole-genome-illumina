#!/usr/bin/env python3

import sys
import argparse
from samplesheet import (
    SampleSheet,
    CORR_SAMPLESHEET,
    DEF_SAMPLE_ID,
    DEF_FW_READS,
    DEF_RV_READS,
    DEF_RUN,
)

RECOGNIZED_COMMANDS = ["read", "write"]


class ReadSampleSheet(object):
    def __init__(self) -> None:
        parser = argparse.ArgumentParser(
            prog="WGS Samplesheet Reader",
            description="Read WGS Samplesheets and update relative paths to samples listed or update previously read samplesheet to absolute paths and save to a masterfile.",
            usage="""process_samplesheet <command> <args> 
            
            The script has two commands:
            - read: Reads samplesheet and looks for each read in the child directories of the samplesheet and updates their relative paths if needed.
            - write: Reads a (processed) samplesheet and updates the paths to absolute paths, then saves the result to a master database file if not already present.
            """,
        )
        parser.add_argument(
            "command",
            help="Subcommand to run. Choice from {}.".format(RECOGNIZED_COMMANDS),
        )
        args = parser.parse_args(sys.argv[1:2])
        if not hasattr(self, args.command):
            print("Unrecognized command")
            parser.print_help()
            exit(1)
        getattr(self, args.command)()

    def read(self):

        parser = argparse.ArgumentParser(
            description="Read samplesheet and correct paths where needed.",
            usage=f"""read_samplesheet.py read -s data/samplesheet.tsv
\ndefaults:
-i: input column header, default: {DEF_SAMPLE_ID}
-f: forward read column header, default:{DEF_FW_READS}
-r: reverse read column header, default: {DEF_RV_READS}
-n: run name, default: {DEF_RUN}
                """,
        )
        parser.add_argument("-s", "--samplesheet", required=True)
        parser.add_argument("-i", "--sample_column", default=DEF_SAMPLE_ID)
        parser.add_argument("-f", "--forward_column", default=DEF_FW_READS)
        parser.add_argument("-r", "--reverse_column", default=DEF_RV_READS)
        parser.add_argument("-n", "--run_name", default=DEF_RUN)
        args = parser.parse_args(sys.argv[2:])

        smpsh = SampleSheet(
            args.samplesheet,
            args.sample_column,
            args.forward_column,
            args.reverse_column,
            None,
            args.run_name,
        )
        smpsh.read_samplesheet()
        smpsh.write_samplesheet()

    def write(self):

        parser = argparse.ArgumentParser(
            description="Read samplesheet, update to absolute paths and store in db.",
            usage="""read_samplesheet.py write -s data/samplesheet.tsv -d database_file.tsv"""
        )
        parser.add_argument("-s", "--samplesheet", default=CORR_SAMPLESHEET)
        parser.add_argument("-d", "--sample_db_dir", required=True)

        args = parser.parse_args(sys.argv[2:])

        smpsh = SampleSheet(
            args.samplesheet, DEF_SAMPLE_ID, DEF_FW_READS, DEF_RV_READS, None, None
        )
        smpsh.read_samplesheet()
        smpsh.update_samplesheet()
        smpsh.update_sampledb()


if __name__ == "__main__":

    ReadSampleSheet()

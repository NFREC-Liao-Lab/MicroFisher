# from ast import arg
from pathlib import Path
import os
import sys


class Config:
    def __init__(self, args, distinct_count=1):
        self.distinct_count = distinct_count
        self.verbose = args.verbose

        try:
            self.min_len = args.min
        except AttributeError:
            self.min_len = 0

        self.centrifuge_path = os.path.expanduser(args.centrifuge_path)

        self.workspace = os.path.expanduser(args.workspace)
        self.out_dir = args.out_dir

        # self.centrifuge_path = args.centrifuge_path
        self.db_path = os.path.join(self.centrifuge_path, args.db_path)
        self.db_name = args.db
        self.db = os.path.join(self.db_path, self.db_name)
        self.threads = args.threads

        self.centrifuge_input_files(args)
        self.centrifuge_output_files()

    def print(self):
        print(self.min_len, self.db_name)

    def centrifuge_input_files(self, args):
        if args.prefix is not None:
            infiles = [os.path.join(
                self.workspace, f"{args.prefix}_R{i}.fastq.gz") for i in [1, 2]]
            self.param_input = f"-1 {infiles[0]} -2 {infiles[1]}"
            self.out_prefix = args.prefix
        elif args.paired is not None:
            infiles = [os.path.join(
                self.workspace, f"{i}") for i in args.paired]
            self.param_input = f"-1 {infiles[0]} -2 {infiles[1]}"
            self.out_prefix = Path(args.paired[0]).stem
        elif args.single is not None:
            infiles = os.path.join(self.workspace, args.single)
            self.param_input = f"-U {infiles}"
            self.out_prefix = Path(args.single).stem

    def centrifuge_output_files(self):
        prefix_outfile = f"result_{self.out_prefix}_min{self.min_len}_db{self.db_name}"
        self.out_output = os.path.join(
            self.out_dir, f"{prefix_outfile}_output.txt")
        self.out_report = os.path.join(
            self.out_dir, f"{prefix_outfile}_report.tsv")
        return self.out_output, self.out_report

    def get_io_files(self):
        return f"{self.param_input} -S {self.out_output} --report-file {self.out_report}"

    def format_params_centrifuge(self):
        params = f"-p {self.threads} -k {self.distinct_count} --min-hitlen {self.min_len} -x {self.db}"
        return params

    def format_params_kreport(self):
        params = f"-x {self.db} {self.out_report}"
        return params

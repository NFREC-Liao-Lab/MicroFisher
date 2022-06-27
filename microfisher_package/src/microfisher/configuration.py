import os


class Config:
    def __init__(self, args, distinct_count=1):
        self.distinct_count = distinct_count
        self.verbose = args.verbose

        try:
            self.min_len = args.min
            self.prefix = args.prefix
        except AttributeError:
            self.min_len = 0
            self.prefix = ""

        self.centrifuge_path = os.path.expanduser(args.centrifuge_path)
        self.workspace = os.path.expanduser(args.workspace)

        self.centrifuge_path = args.centrifuge_path
        self.db_path = os.path.join(self.centrifuge_path, args.db_path)
        self.db_name = args.db
        self.db = os.path.join(self.db_path, self.db_name)
        self.threads = args.threads

    def print(self):
        print(self.min_len, self.db_name)

    def centrifuge_output_files(self):
        prefix_outfile = f"result_{self.prefix}_min{self.min_len}_db{self.db_name}"
        self.out_output = os.path.join(self.workspace, f"{prefix_outfile}_output.txt")
        self.out_report = os.path.join(self.workspace, f"{prefix_outfile}_report.tsv")
        return self.out_output, self.out_report

    def format_params_centrifuge(self):
        params = f"-p {self.threads} -k {self.distinct_count} --min-hitlen {self.min_len} -x {self.db}"
        return params

    def format_params_kreport(self):
        params = f"-x {self.db} {self.out_report}"
        return params

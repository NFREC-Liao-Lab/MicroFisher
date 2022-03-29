import os
import sys
import subprocess
import shlex
from config import Config

# del sys.modules["config"]
# import importlib
# importlib.reload(config.Config)

# cc = Config(3,min_len=10)
# cc.print()
# cc.cpus
# cc.distinct_count
# cc.min_len


class Centrifuge:

    def __init__(self, config):
        self.config = config
        # self.centrifuge_path = os.path.expanduser(args.centrifuge_path)
        # self.workspace = os.path.expanduser(args.workspace)
        # self.prefix = args.prefix# f"simulating_{num}species_r{replicate}"
        #
        # self.config = Config(cpus=8, distinct_count=1, min_len=args.min,
        #             centrifuge_path=self.centrifuge_path
        #             db_path=args.db_path, db_name=args.db
        #             )

        infile = [os.path.join(
            self.config.workspace, f"{self.config.prefix}_R{i}.fastq.gz") for i in [1, 2]]
        out_result, out_result = config.centrifuge_output_files()
        params = self.config.format_params_centrifuge()
        config_files = f"-1 {infile[0]} -2 {infile[1]} -S {out_result} --report-file {self.config.out_report}"
        prog = os.path.join(self.config.centrifuge_path, "centrifuge")
        commands = f"{prog} {params} {config_files}"
        self.command_list = shlex.split(commands)
        if config.verbose > 0:
            print(f"==DEBUG== Execute commands: {commands}")

    def run(self):
        try:
            result = subprocess.run(self.command_list)
        except Exception as err:
            print(f"Unexpected {err=}, {type(err)=}")

        is_output_exist = os.path.exists(self.config.out_report)
        if result.returncode == 0 and is_output_exist:
            print(f"==DEBUG== Centrifuge complete.")
            return True
        return False

    # prog_k = "centrifuge-kreport"
    # out_kreport = f"{workspace}/result_{prefix_out}_kreport.tsv"
    # config = f"-x {cc.db_name} {out_report} {out_kreport}"

    # command_list = shlex.split("ls -l -F --color='always'")
    # subprocess.run(command_list)

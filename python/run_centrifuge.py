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


    def __init__(self, args):
        num=100
        replicate=5
        self.workspace = os.path.expanduser(args.workspace)
        centrifuge_path = os.path.expanduser(args.centrifuge)
        db_name = args.db
        db_path = os.path.join(centrifuge_path, args.db_path)
        cc = Config(cpus=4, distinct_count=1, min_len=120,
                    db_path=db_path, db_name=db_name
                    )


        prefix = f"simulating_{num}species_r{replicate}"
        prefix_infile = f"{prefix}.short_read"
        infile = [os.path.join(self.workspace, f"{prefix_infile}_R{i}.fastq.gz") for i in [1,2]]

        prefix_outfile = f"result_{prefix}_{cc.min_len}_{cc.db_name}"
        out_result = os.path.join(self.workspace, f"{prefix_outfile}_output.txt")
        out_report = os.path.join(self.workspace, f"{prefix_outfile}_report.tsv")

        params = cc.format()
        config_files = f"-1 {infile[0]} -2 {infile[1]} -S {out_result} --report-file {out_report}"
        prog = "centrifuge"
        commands = f"{prog} {params} {config_files}"
        command_list = shlex.split(commands)
        if (args.verbose > 0):
            print(f"Execute commands: {commands}")
        subprocess.run(command_list)


    # prog_k = "centrifuge-kreport"
    # out_kreport = f"{workspace}/result_{prefix_out}_kreport.tsv"
    # config = f"-x {cc.db_name} {out_report} {out_kreport}"


    # command_list = shlex.split("ls -l -F --color='always'")
    # subprocess.run(command_list)

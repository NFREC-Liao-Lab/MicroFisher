import os
import sys
import subprocess
import shlex
# from .config import Config

# del sys.modules["config"]
# import importlib
# importlib.reload(config.Config)

# cc = Config(3,min_len=10)
# cc.print()
# cc.cpus
# cc.distinct_count
# cc.min_len


class CentrifugeKReport:

    def __init__(self, config):
        self.config = config

        prefix, ext = os.path.splitext(config.out_report)
        self.out_kreport = f"{prefix}_kreport.tsv"

        prog = os.path.join(config.centrifuge_path, "centrifuge-kreport")
        params = config.format_params_kreport()
        commands = f"{prog} {params}"
        self.command_list = shlex.split(commands)
        if config.verbose > 0:
            print(f"==DEBUG== Execute commands: {commands}")

    def run(self):
        try:
            kreport_output = subprocess.run(self.command_list,
                                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)  # capture_output=True)
            if kreport_output.stdout:
                with open(self.out_kreport, "wb") as f:
                    f.write(kreport_output.stdout)
            if self.config.verbose > 1 and kreport_output.stderr:
                with open(self.out_kreport + ".err", "wb") as f:
                    f.write(kreport_output.stderr)

        except Exception as err:
            print(f"Unexpected {err}, {type(err)}")
        # print(kreport_output)
        # print(kreport_output.stdout)


# /Users/steven/workspace/centrifuge/centrifuge -x /Users/steven/workspace/centrifuge/example/index/test -1 /Users/steven/workspace/evol1_R1.fastq.gz -2 /Users/steven/workspace/evol1_R2.fastq.gz -S /Users/steven/workspace/result_evol1_min120_dbtest_output.txt --report-file /Users/steven/workspace/result_evol1_min 120_dbtest_report.tsv
    # prog_k = "centrifuge-kreport"
    # out_kreport = f"{workspace}/result_{prefix_out}_kreport.tsv"
    # config = f"-x {cc.db_name} {out_report} {out_kreport}"

    # command_list = shlex.split("ls -l -F --color='always'")
    # subprocess.run(command_list)

import os
import subprocess
import shlex
# from .config import Config


class CentrifugeKReport:

    def __init__(self, config):
        self.config = config

        prefix, ext = os.path.splitext(config.out_report)
        self.out_kreport = f"{prefix}_kreport.tsv"

        prog = os.path.join(config.centrifuge_path, "centrifuge-kreport")
        params = config.format_params_kreport()
        commands = f"{prog} {params}"
        self.command_list = shlex.split(commands)
        if config.verbose > 1:
            print(f"===DEBUG=== Execute commands: {commands}")

    def run(self):
        try:
            kreport_output = subprocess.run(self.command_list, stdout=subprocess.PIPE,
                                            stderr=subprocess.PIPE)  # capture_output=True)
            if kreport_output.stdout:
                with open(self.out_kreport, "wb") as f:
                    f.write(kreport_output.stdout)
            if self.config.verbose > 1 and kreport_output.stderr:
                with open(self.out_kreport + ".err", "wb") as f:
                    f.write(kreport_output.stderr)
            return True
        except Exception as err:
            print(f"Unexpected {err}, {type(err)}")
        # print(kreport_output)
        # print(kreport_output.stdout)

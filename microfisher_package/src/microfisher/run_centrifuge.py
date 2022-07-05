import os
import subprocess
import shlex
# from .config import Config

# del sys.modules["config"]
# import importlib
# importlib.reload(config.Config)


class Centrifuge:

    def __init__(self, config):
        self.config = config

        params = self.config.format_params_centrifuge()
        config_files = self.config.get_io_files()

        prog = os.path.join(self.config.centrifuge_path, "centrifuge")
        commands = f"{prog} {params} {config_files}"
        self.command_list = shlex.split(commands)
        if self.config.verbose > 0:
            print(f"==DEBUG== Centrifuge: Start.")
        if self.config.verbose > 1:
            print(f"===DEBUG=== Execute commands: {commands}")

    def run(self):
        try:
            result = subprocess.run(self.command_list)
        except Exception as err:
            print(f"Unexpected {err}, {type(err)}")

        is_output_exist = os.path.exists(self.config.out_report)
        if is_output_exist and result.returncode == 0:
            print("==DEBUG== Centrifuge: Done.")
            return True
        return False

    # prog_k = "centrifuge-kreport"
    # out_kreport = f"{workspace}/result_{prefix_out}_kreport.tsv"
    # config = f"-x {cc.db_name} {out_report} {out_kreport}"

    # command_list = shlex.split("ls -l -F --color='always'")
    # subprocess.run(command_list)

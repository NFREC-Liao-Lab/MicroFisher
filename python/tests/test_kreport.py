import pytest
import shlex
from tests.test_config import setup_args, config
# from config import Config
from run_kreport import CentrifugeKReport


@pytest.mark.usefixtures("setup_args", "config")
def test_kreport(config):
    config.centrifuge_output_files()
    kreport = CentrifugeKReport(config)
    commands = ("cpath/centrifuge-kreport "
    "-x cpath/db/dbName workspace/result_example_min100_dbdbName_report.tsv"
                )
    expected = shlex.split(commands)
    assert kreport.command_list == expected

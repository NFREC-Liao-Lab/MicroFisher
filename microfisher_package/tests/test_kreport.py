import shlex

import pytest

from microfisher.run_kreport import CentrifugeKReport
from tests.test_config import config, setup_args


@pytest.mark.usefixtures("setup_args", "config")
def test_kreport(config):
    config.centrifuge_output_files()
    kreport = CentrifugeKReport(config)
    commands = ("cpath/centrifuge-kreport "
                "-x cpath/db/dbName workspace/result_example_min100_dbdbName_report.tsv"
                )
    expected = shlex.split(commands)
    assert kreport.command_list == expected

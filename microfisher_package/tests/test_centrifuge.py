import pytest
import shlex
from tests.test_config import setup_args, config
from src.microfisher.run_centrifuge import Centrifuge


@pytest.mark.usefixtures("setup_args", "config")
def test_centrifuge(config):
    centrifuge = Centrifuge(config)
    commands = ("cpath/centrifuge -p 4 -k 1 --min-hitlen 100 -x "
                "cpath/db/dbName "
                "-1 workspace/example_R1.fastq.gz -2 workspace/example_R2.fastq.gz "
                "-S workspace/merged_results/result_example_min100_dbdbName_output.txt "
                "--report-file workspace/merged_results/result_example_min100_dbdbName_report.tsv "
                )
    expected = shlex.split(commands)
    assert centrifuge.command_list == expected

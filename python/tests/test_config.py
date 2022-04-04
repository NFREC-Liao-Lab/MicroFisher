import argparse
import pytest
from config import Config
import tempfile


@pytest.fixture
def setup_args():
    parser = argparse.ArgumentParser()
    parser.set_defaults(verbose=1, min=100, prefix="example",
        workspace="workspace",
        centrifuge_path="cpath", db_path="db", db="dbName")
    args = parser.parse_args()
    return args


@pytest.fixture
def config(setup_args):
    return Config(setup_args)


def test_init(setup_args):
    config = Config(setup_args)
    assert config.verbose == 11
    assert config.min_len == 1001
    assert config.prefix == "example"
    assert config.centrifuge_path == "cpath"
    assert config.workspace == "workspace"
    assert config.db == "cpath/db/dbName"


def test_output_files(config):
    output = config.centrifuge_output_files()
    expected = ("workspace/result_example_min100_dbdbName_output.txt", "workspace/result_example_min100_dbdbName_report.tsv")
    assert expected == output


def test_format(config):
    config.centrifuge_output_files()
    output = config.format_params_centrifuge()
    assert "-p 1 -k 1 --min-hitlen 100 -x cpath/db/dbName" == output

    output = config.format_params_kreport()
    assert "-x cpath/db/dbName workspace/result_example_min100_dbdbName_report.tsv" == output

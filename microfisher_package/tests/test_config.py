import argparse
import pytest
import tempfile
from src.microfisher.configuration import Config
from src.microfisher import microfisher


@pytest.fixture
def setup_args():
    parser = argparse.ArgumentParser()
    parser.set_defaults(verbose=1, min=100, prefix="example",
                        threads=4,
                        workspace="workspace",
                        out_dir="merged_results",
                        centrifuge_path="cpath", db_path="db", db="dbName")
    args = parser.parse_args([])
    args = microfisher.parse_output_dir(args)
    return args


@pytest.fixture
def config(setup_args):
    return Config(setup_args)


def test_init(setup_args):
    config = Config(setup_args)
    assert config.verbose == 1
    assert config.min_len == 100
    assert config.threads == 4
    assert config.centrifuge_path == "cpath"
    assert config.workspace == "workspace"
    assert config.db == "cpath/db/dbName"
    assert config.param_input == "-1 workspace/example_R1.fastq.gz -2 workspace/example_R2.fastq.gz"
    assert config.out_prefix == "example"


def test_centrifuge_input_single(setup_args):
    setup_args.prefix = None
    setup_args.paired = None
    setup_args.single = "single_end.fastq.gz"
    config = Config(setup_args)
    assert config.workspace == "workspace"
    assert config.db == "cpath/db/dbName"
    assert config.param_input == "-U workspace/single_end.fastq.gz"
    assert config.out_prefix == "single_end.fastq"


def test_centrifuge_input_paired(setup_args):
    setup_args.prefix = None
    setup_args.paired = ["paired_1.fastq.gz", "paired_2.fastq.gz"]
    setup_args.single = None
    config = Config(setup_args)
    assert config.workspace == "workspace"
    assert config.db == "cpath/db/dbName"
    assert config.param_input == "-1 workspace/paired_1.fastq.gz -2 workspace/paired_2.fastq.gz"
    assert config.out_prefix == "paired_1.fastq"


def test_output_files(config):
    output = config.centrifuge_output_files()
    expected = ("workspace/merged_results/result_example_min100_dbdbName_output.txt",
                "workspace/merged_results/result_example_min100_dbdbName_report.tsv")
    assert expected == output


def test_format(config):
    config.centrifuge_output_files()
    output = config.format_params_centrifuge()
    assert "-p 4 -k 1 --min-hitlen 100 -x cpath/db/dbName" == output

    output = config.format_params_kreport()
    assert "-x cpath/db/dbName workspace/merged_results/result_example_min100_dbdbName_report.tsv" == output

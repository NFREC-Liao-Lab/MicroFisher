import argparse
import pytest
from config import Config

@pytest.fixture
def config_args():
    parser = argparse.ArgumentParser()
    parser.set_defaults(verbose=1, min=100, prefix="example",
        workspace="workspace",
        centrifuge_path="cpath", db_path="db", db="dbName")
    args = parser.parse_args()
    return args

def test_init(config_args):

    config = Config(config_args)
    assert config.verbose == 1
    assert config.min_len == 100
    assert config.prefix == "example"
    assert config.centrifuge_path == "cpath"
    assert config.workspace == "workspace"
    assert config.db == "cpath/db/dbName"

def test_init(config_args):

    config = Config(config_args)
    assert config.verbose == 1
    assert config.min_len == 100
    assert config.prefix == "example"
    assert config.centrifuge_path == "cpath"
    assert config.workspace == "workspace"
    assert config.db == "cpath/db/dbName"

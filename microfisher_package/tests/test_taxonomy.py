import pytest
import tempfile
from src.microfisher import taxonomy
from src.microfisher import taxonomy_utils

# @pytest.fixture
# def setup_args():
#     parser = argparse.ArgumentParser()
#     parser.set_defaults(verbose=1, min=100, prefix="example",
#         workspace="workspace",
#         centrifuge_path="cpath", db_path="db", db="dbName")
#     args = parser.parse_args()
#     return args
#
#
# @pytest.fixture
# def config(setup_args):
#     return Config(setup_args)

def test_taxonomy_desired():
    taxid = 9606
    actual = taxonomy.get_desired_taxa_ranks(taxid, desired_ranks=taxonomy.DESIRED_RANKS)
    expected = {'family': (9604, 'Hominidae'), 'genus': (9605, 'Homo'), 'species': (9606, 'Homo sapiens')}
    assert actual == expected


def test_taxonomy_default():
    taxid = 9606
    actual = taxonomy.get_desired_taxa_ranks(taxid)# =["phylum", "class", "order", "family", "genus", "species"])
    expected = {'kingdom': (33208, 'Metazoa'), 'phylum': (7711, 'Chordata'), 'class': (40674, 'Mammalia'), 'order': (9443, 'Primates'), 'family': (
        9604, 'Hominidae'), 'genus': (9605, 'Homo'), 'species': (9606, 'Homo sapiens')}
    assert actual == expected


def test_taxonomy_1():
    taxid = 9606
    actual = taxonomy.get_desired_taxa_ranks(taxid, desired_ranks=["species"])
    expected = {'species': (9606, 'Homo sapiens')}
    assert actual == expected

    actual = taxonomy.get_desired_taxa_ranks(taxid, desired_ranks=["family", "species"])
    expected = {'family': (9604, 'Hominidae'), 'species': (9606, 'Homo sapiens')}
    assert actual == expected


def test_summarise_data():
    data = {"D1": [
        {"genus": "genus_1", "species": "species_1", "numReads": 10},
        {"genus": "genus_1", "species": "species_1", "numReads": 11},
        {"genus": "genus_1", "species": "species_2", "numReads": 12},
        {"genus": "genus_2", "species": "species_3", "numReads": 13},
        {"genus": "genus_2", "species": "species_4", "numReads": 14}
    ],
        "D2": [
        {"genus": "genus_1", "species": "species_1", "numReads": 30},
        {"genus": "genus_1", "species": "species_1", "numReads": 31},
        {"genus": "genus_2", "species": "species_2", "numReads": 32},
        {"genus": "genus_2", "species": "species_3", "numReads": 33},
        {"genus": "genus_1", "species": "species_4", "numReads": 34}
    ]
    }
    actual = taxonomy_utils.summarise_data(data["D1"], rank="species")
    expected = {"species_1": 21, "species_2": 12,
                "species_3": 13, "species_4": 14}
    assert actual == expected

    actual = taxonomy_utils.summarise_data(data["D1"], rank="genus")
    expected = {"genus_2": 27, "genus_1": 33}
    assert actual == expected

    actual = taxonomy_utils.summarise_data(data["D2"], rank="genus")
    expected = {"genus_2": 65, "genus_1": 95}
    assert actual == expected


def test_normalise_data_1():

    data = {"a": 13, "b": 31, "c": 29, "d": 27}
    actual = taxonomy_utils.normalise_data(data)
    expected = {"a": 0.13, "b": 0.31, "c": 0.29, "d": 0.27}
    assert actual == expected


def test_normalise_data_2():

    data = {"a": 260, "b": 620, "c": 580, "d": 540}
    actual = taxonomy_utils.normalise_data(data)
    expected = {"a": 0.13, "b": 0.31, "c": 0.29, "d": 0.27}
    assert actual == expected


def test_inverse_normalise_data_2():

    data = {"a": 13, "b": 31, "c": 29, "d": 27}
    actual = taxonomy_utils.inverse_normalise_data(data)
    expected = {"a": 0.4256927, "b": 0.1785163, "c": 0.1908278, "d": 0.2049632}
    assert actual == pytest.approx(expected)


def test_weight_data():

    data = {"a": 10, "b": 20, "c": 30.30, "d": -40}
    weight = 0.5
    actual = taxonomy_utils.weight_data(data, weight)
    expected = {"a": 5, "b": 10, "c": 15.15, "d": -20}
    assert actual == expected

    weight = 2.2
    actual = taxonomy_utils.weight_data(data, weight)
    expected = {"a": 22, "b": 44, "c": 66.66, "d": -88}
    assert actual == pytest.approx(expected)


def test_create_length_dict():
    file_name = ["F1", "F2", "F3"]
    actual = taxonomy_utils.create_length_dict(file_name, None)
    assert actual is None
    length = [12, 10009, -3.2]
    actual = taxonomy_utils.create_length_dict(file_name, length)
    print(actual)
    expected = {"F1": 12, "F2": 10009, "F3": -3.2}
    assert actual == expected

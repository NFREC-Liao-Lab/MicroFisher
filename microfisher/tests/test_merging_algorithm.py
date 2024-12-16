import pytest
import src.microfisher.merging_algorithm as malg
from src.microfisher import taxonomy_utils

data_summary = {
    "data1": {
        "s1": 11, "s2": 21, "s3": 31, "s4": 41
    },
    "data2": {
        "s1": 12, "s2": 22, "s3": 32, "s5": 52
    },
    "data3": {
        "s1": 13, "s2": 23, "s4": 43, "s6": 63
    }
}


def test_merging_boolean():

    actual = malg.a_boolean(data_summary, 1)
    expected = {"s1": 3, "s2": 3, "s3": 2, "s4": 2, "s5": 1, "s6": 1}
    assert actual == expected

    actual = malg.a_boolean(data_summary, 2)
    expected = {"s1": 3, "s2": 3, "s3": 2, "s4": 2}
    assert actual == expected

    actual = malg.a_boolean(data_summary, 3)
    expected = {"s1": 3, "s2": 3}
    assert actual == expected


def test_merging_equal():

    actual = malg.a_equal(data_summary)
    expected = {"s1": 36, "s2": 66, "s3": 63, "s4": 84, "s5": 52, "s6": 63}
    assert actual == expected


def test_cal_percentage_1():

    data = {"s1": 36, "s2": 66, "s3": 63, "s4": 84, "s5": 52, "s6": 63}
    actual = malg.calculate_percentage(data)
    expected = {k: v/364 for k, v in data.items()}
    assert actual == expected


def test_cal_percentage_sq():

    data = {i: i*i for i in range(100)}
    total = sum(data.values())
    expected = {k: v/total for k, v in data.items()}
    actual = malg.calculate_percentage(data)
    assert actual == expected


data_summary_2 = {
    "data1": {
        "s1": 10, "s2": 20, "s3": 70
    },
    "data2": {
        "s1": 10, "s3": 30, "s4": 60
    },
    "data3": {
        "s1": 10, "s2": 50, "s4": 40
    }
}


def test_cal_percentage_weighted_basic():

    actual = malg.a_weighted(data_summary_2, cent_length_dict=None,
                             db_mean_length=None)
    expected = {"s1": 0.1, "s2": 0.7/3, "s3": 1/3, "s4": 1/3}
    assert actual == pytest.approx(expected)


def test_cal_percentage_weighted_cent():

    cent_length = {
        "data1": 100,
        "data2": 50,
        "data3": 200,
    }

    actual = malg.a_weighted(data_summary_2, cent_length_dict=cent_length,
                             db_mean_length=None)
    expected = {"s1": 0.1, "s2": 0.2*2/7 + 0.5*4/7,
                "s3": 0.7*2/7 + 0.3/7,
                "s4": 0.6/7+0.4*4/7}
    assert actual == pytest.approx(expected)


def test_cal_percentage_weighted_db():

    db_length = {
        "data1": 100,
        "data2": 50,
        "data3": 200,
    }
    actual = malg.a_weighted(data_summary_2, cent_length_dict=None,
                             db_mean_length=db_length)
    inv_total = 1/100 + 1/50 + 1/200
    inv_prop = [(1/100)/inv_total, (1/50)/inv_total, (1/200)/inv_total]
    expected = {"s1": 0.1, "s2": 0.2*inv_prop[0] + 0.5*inv_prop[2],
                "s3": 0.7*inv_prop[0] + 0.3*inv_prop[1],
                "s4": 0.6*inv_prop[1] + 0.4*inv_prop[2]}
    assert actual == pytest.approx(expected)


def test_cal_percentage_weighted_db2():

    cent_length = {
        "data1": 100,
        "data2": 200,
        "data3": 300,
    }
    db_length = {
        "data1": 100,
        "data2": 50,
        "data3": 200,
    }
    actual = malg.a_weighted(data_summary_2, cent_length_dict=cent_length,
                             db_mean_length=db_length)
    inv_total = 1/100 + 1/50 + 1/200
    inv_prop = [(1/100)/inv_total, (1/50)/inv_total, (1/200)/inv_total]
    prop = [inv_prop[0]*1/6, inv_prop[1]*2/6, inv_prop[2]*3/6]
    expected_temp = {
        "s1": sum([0.1 * p for p in prop]),
        "s2": 0.2 * prop[0] + 0.5 * prop[2],
        "s3": 0.7 * prop[0] + 0.3 * prop[1],
        "s4": 0.6 * prop[1] + 0.4 * prop[2],
    }
    expected = taxonomy_utils.normalise_data(expected_temp)
    assert actual == pytest.approx(expected)

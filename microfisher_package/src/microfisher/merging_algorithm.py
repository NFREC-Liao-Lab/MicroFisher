# Merging algorithms
from . import taxonomy_utils

MODE_CHOICES = ["weighted", "boolean", "raw", "weighted_abundance_only",
                "weighted_centlength_only"]
DB_LIST = {"ITS+LSU": ["ITS1_fisher", "ITS2_fisher", "LSU_D1_fisher_new", "LSU_D2_fisher_new"],
           "ITS": ["ITS1_fisher", "ITS2_fisher"],
           "LSU": ["LSU_D1_fisher_new", "LSU_D2_fisher_new"]}
FILTER_DEFAULT = 0.00001
DEFAULT_DB_LENGTH = {
    "ITS1": 188.8,
    "ITS2": 186.8,
    "LSUD1": 172.8,
    "LSUD2": 189.4
}


def a_boolean(data_summary, min_db_count):

    counts = dict()
    for data_each in data_summary.values():
        for k, v in data_each.items():
            try:
                counts[k] += 1
            except KeyError:
                counts[k] = 1

    results = {k: v for k, v in counts.items() if v >= min_db_count}
    return(results)


def a_equal(data_summary):

    results = dict()
    for data_each in data_summary.values():
        for k, v in data_each.items():
            try:
                results[k] += v
            except KeyError:
                results[k] = v
    return(results)


def a_weighted(data_summary, cent_length_dict=None, db_mean_length=None):
    # print(cent_length_dict)
    # print(db_mean_length)
    data_mode = {k: taxonomy_utils.normalise_data(v)
                 for k, v in data_summary.items()}

    if cent_length_dict is not None:
        cent_weight = taxonomy_utils.normalise_data(cent_length_dict)
        print("Cent:", cent_weight)
        data_mode = {k: taxonomy_utils.weight_data(v, cent_weight[k])
                     for k, v in data_mode.items()}

    if db_mean_length is not None:
        db_weight = taxonomy_utils.inverse_normalise_data(db_mean_length)
        print("DB:", db_weight)
        data_mode = {k: taxonomy_utils.weight_data(v, db_weight[k])
                     for k, v in data_mode.items()}

    results = a_equal(data_mode)
    # normalised_factor = sum(results.values())
    # percentage = {k: v / normalised_factor for k, v in results.items()}
    percentage = taxonomy_utils.normalise_data(results)
    return(percentage)


def a_probability(data_summary, report_length_dict):

    # data_mode = {k:taxonomy_utils.normalise_data(v) for k, v in data_summary.items()}
    # # if report_length_dict is not None:
    # #     max_length = max(report_length_dict.values())
    # #     weight_length = {k:v/max_length for k, v in report_length_dict.items()}
    # #     data_mode = {k:taxonomy_utils.weight_data(v, weight_length[k]) for k, v in data_mode.items()}
    # results = a_raw(data_mode)
    # normalised_factor = len(data_summary)
    # percentage = {k:v/normalised_factor for k, v in results.items()}
    probability = None
    return(probability)


def calculate_percentage(results):
    total = sum(results.values())
    percentage = {k: v / total for k, v in results.items()}
    return(percentage)

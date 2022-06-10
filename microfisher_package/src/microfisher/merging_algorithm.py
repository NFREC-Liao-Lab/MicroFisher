# Merging algorithms
from . import taxonomy_utils

MODE_CHOICES = ["boolean", "raw", "weighted", "weighted_length", "probability"]


def a_boolean(data_summary, min_db_conut):

    counts = dict()
    for data_each in data_summary.values():
        for k, v in data_each.items():
            try:
                counts[k] += 1
            except KeyError:
                counts[k] = 1

    results = {k: v for k, v in counts.items() if v >= min_db_conut}
    return(results)


def a_raw(data_summary):

    results = dict()
    for data_each in data_summary.values():
        for k, v in data_each.items():
            try:
                results[k] += v
            except KeyError:
                results[k] = v
    return(results)


def a_weighted(data_summary, report_length_dict=None):

    data_mode = {k: taxonomy_utils.normalise_data(v)
                 for k, v in data_summary.items()}
    if report_length_dict is not None:
        max_length = max(report_length_dict.values())
        weight_length = {k: v / max_length for k,
                         v in report_length_dict.items()}
        data_mode = {k: taxonomy_utils.weight_data(v, weight_length[k])
                     for k, v in data_mode.items()}
    results = a_raw(data_mode)
    normalised_factor = len(data_summary)
    percentage = {k: v / normalised_factor for k, v in results.items()}
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

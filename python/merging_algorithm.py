# Merging algorithms
import taxanomy_utils

MODE_CHOICES = ["raw", "weighted", "probability"]


def a_raw(data_mode):

    results = dict()
    for data_each in data_mode.values():
        for k, v in data_each.items():
            try:
                results[k] += v
            except KeyError as e:
                results[k] = v
    return(results)



def a_weighted(data_summary, report_length_dict):

    data_mode = {k:taxanomy_utils.normalise_data(v) for k, v in data_summary.items()}
    if report_length_dict is not None:
        max_length = max(report_length_dict.values())
        weight_length = {k:v/max_length for k, v in report_length_dict.items()}
        data_mode = {k:taxanomy_utils.weight_data(v, weight_length[k]) for k, v in data_mode.items()}
    results = a_raw(data_mode)
    normalised_factor = len(data_summary)
    percentage = {k:v/normalised_factor for k, v in results.items()}
    return(percentage)


def calculate_percentage(results):
    total = sum(results.values())
    percentage = {k:v/total for k, v in results.items()}
    return(percentage)

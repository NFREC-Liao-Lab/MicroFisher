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



def a_weighted(data_summary):

    data_mode = {k:taxanomy_utils.normalise_data(v) for k, v in data_summary.items()}
    results = a_raw(data_mode)
    normalised_factor = len(data_summary)
    percentage = {k:v/normalised_factor for k, v in results.items()}
    return(percentage)


def calculate_percentage(results):
    total = sum(results.values())
    percentage = {k:v/total for k, v in results.items()}
    return(percentage)

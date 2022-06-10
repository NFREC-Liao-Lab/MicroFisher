
def format_merge_single(results, threshold, label="proportion"):

    sorted_key = sorted(results, key=results.get, reverse=True)
    output_header = f"taxID\tname\t{label}\n"

    output_list = [f"{k[0]}\t{k[1]}\t{results[k]:g}" for k in sorted_key]
    output_list = output_header + "\n".join(output_list)

    output_filter_list = None
    if threshold is not None:
        filtered_key = [k for k in sorted_key if results[k] > threshold]
        output_filter_list = [f"{k[0]}\t{k[1]}\t{results[k]:g}" for k in filtered_key]
        output_filter_list = output_header + "\n".join(output_filter_list)

    return(output_list, output_filter_list)


def format_merge_two(results, percentage, threshold, label=["merged_reads", "percentage"]):

    sorted_key = sorted(results, key=results.get, reverse=True)

    output_header = f"taxID\tname\t{label[0]}\t{label[1]}\n"
    output_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}" for k in sorted_key]
    output_list = output_header + "\n".join(output_list)

    filtered_key = [k for k in sorted_key if percentage[k] > threshold]
    output_filter_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}"
                          for k in filtered_key]
    output_filter_list = output_header + "\n".join(output_filter_list)

    return(output_list, output_filter_list)

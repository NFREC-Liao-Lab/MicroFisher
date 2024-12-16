from . import taxonomy


def generate_taxa_full_rank_str(taxid):
    full_ranks = taxonomy.get_desired_taxa_ranks(taxid)
    rank_list = [full_ranks[t][1]
                 if t in full_ranks else "None" for t in taxonomy.FULL_RANKS]
    rank_str = "\t".join(rank_list)
    return rank_str


def generate_single_results_str(key, results):
    rank_str = generate_taxa_full_rank_str(key[0])
    output = f"{key[0]}\t{key[1]}\t{results[key]:g}\t{rank_str}\n"
    return output


def generate_merge_two_results_str(key, results, percentage):
    rank_str = generate_taxa_full_rank_str(key[0])
    output = f"{key[0]}\t{key[1]}\t{results[key]}\t{percentage[key]:.6f}\t{rank_str}\n"
    return output


def format_merge_single(results, threshold, label="proportion"):
    sorted_key = sorted(results, key=results.get, reverse=True)

    output_header = f"taxID\tname\t{label}\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
    all_output_dict = {
        key: generate_single_results_str(
            key, results) for key in sorted_key}
    output_list = [output_header]
    output_list.extend(list(all_output_dict.values()))

    output_filter_list = [output_header]
    if threshold is not None:
        filtered_key = [k for k in sorted_key if results[k] > threshold]
        filtered_output_list = [all_output_dict[k] for k in filtered_key]
        output_filter_list.extend(filtered_output_list)

    return (output_list, output_filter_list)


def format_merge_two(results, percentage, threshold,
                     label=["merged_reads", "percentage"]):

    sorted_key = sorted(results, key=results.get, reverse=True)

    output_header = f"taxID\tname\t{label[0]}\t{label[1]}\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
    all_output_dict = {key: generate_merge_two_results_str(
        key, results, percentage) for key in sorted_key}
    output_list = [output_header]
    output_list.extend(list(all_output_dict.values()))

    output_filter_list = [output_header]
    filtered_key = [k for k in sorted_key if percentage[k] > threshold]
    filtered_output_list = [all_output_dict[k] for k in filtered_key]
    output_filter_list.extend(filtered_output_list)

    return (output_list, output_filter_list)

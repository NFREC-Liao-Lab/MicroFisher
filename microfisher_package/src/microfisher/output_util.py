import sys

from src.microfisher import taxonomy


def generate_taxa_full_rank_str(taxid):
    
    full_ranks = taxonomy.get_desired_taxa_ranks(taxid)
    rank_list = [full_ranks[t][1] if t in full_ranks else "None" for t in taxonomy.FULL_RANKS ]
    rank_str = "\t".join(rank_list)
    return rank_str

def generate_single_results(key, results):

    rank_str = generate_taxa_full_rank_str(key[0])
    output = f"{key[0]}\t{key[1]}\t{results[key]:g}\t{rank_str}"
    return output

def format_merge_single(results, threshold, label="proportion"):

    sorted_key = sorted(results, key=results.get, reverse=True)
    # output_header = f"taxID\tname\t{label}\n"
    output_header = f"taxID\tname\t{label}\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\n"
    all_output_dict = {key: generate_single_results(key, results) for key in sorted_key}
    
    output_list = [output_header]

    output_list.extend(list(all_output_dict.values()))

    output_filter_list = [output_header]
    if threshold is not None:
        filtered_key = [k for k in sorted_key if results[k] > threshold]
        # print(filtered_key, threshold)
        # print(all_output[filtered_key])
        # output_filter_list = [f"{k[0]}\t{k[1]}\t{results[k]:g}" for k in filtered_key]
        filtered_output_list  = [all_output_dict[k] for k in filtered_key]
        
        output_filter_list.extend(filtered_output_list)

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

import os
import sys

from . import taxonomy_utils
from . import merging_algorithm
from . import output_util


def core_merging_data(mode, parsed_data, rank, threshold, cent_length_dict=None, db_length_dict=None):

    try:
        assert mode in merging_algorithm.MODE_CHOICES
    except AssertionError as e:
        sys.exit(f"Invalided combining mode: {mode}. {e}")

    data_summary = {k: taxonomy_utils.summarise_data(v, rank)
                    for k, v in parsed_data.items()}

    if mode == "raw":
        results = merging_algorithm.a_equal(data_summary)
        percentage = merging_algorithm.calculate_percentage(results)
        output_list, output_filter_list = output_util.format_merge_two(
            results, percentage, threshold=threshold, label=["merged_reads", "percentage"])

    elif mode == "boolean":
        results = merging_algorithm.a_boolean(data_summary, min_db_count=2)
        # sorted_key = sorted(results, key=results.get, reverse=True)
        output_list, output_filter_list = output_util.format_merge_single(
            results, threshold=None, label="db_conut")


    elif mode == "weighted":
        results = merging_algorithm.a_weighted(
            data_summary, cent_length_dict, db_length_dict)
        output_list, output_filter_list = output_util.format_merge_single(
            results, threshold=threshold, label="proportion")

    elif mode == "weighted_abundance_only":
        results = merging_algorithm.a_weighted(data_summary, None, None)
        output_list, output_filter_list = output_util.format_merge_single(
            results, threshold=threshold, label="proportion")

    elif mode == "weighted_centlength_only":
        results = merging_algorithm.a_weighted(
            data_summary, cent_length_dict, None)
        output_list, output_filter_list = output_util.format_merge_single(
            results, threshold=threshold, label="proportion")

    elif mode == "probability":
        sys.exit("Exit: mode not yet implemented")
        results = merging_algorithm.a_probability(data_summary)

    return output_list, output_filter_list


def run(args):
    report_files = [os.path.join(args.workspace, f) for f in args.combine]
    out_dir = os.path.join(args.workspace, args.out_dir)
    out_prefix = args.out_prefix
    threshold = args.filter
    mode = args.mode
    desired_ranks = args.desired_ranks
    # min_db_conut = args.min_overlap
    # length_list = args.length

    parsed_data = taxonomy_utils.parse_report_files(report_files)
    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass

    cent_length_dict = taxonomy_utils.create_length_dict(
        report_files, args.cent_length)
    db_length_dict = taxonomy_utils.create_length_dict(
        report_files, args.db_length)
    # print(db_length_dict)
    if args.verbose > 0:
        print("==DEBUG== Combining reports: START.")
    for rank in desired_ranks:
        if args.verbose > 1:
            print(f"===DEBUG=== Combining reports for {rank}.")
        output_list, output_filter_list = core_merging_data(
            mode, parsed_data, rank, threshold, cent_length_dict, db_length_dict)

        outfile = os.path.join(out_dir, f"{out_prefix}_taxa_{rank}.tsv")
        with open(outfile, "w") as fout:
            fout.writelines(output_list)

        outfile = os.path.join(
            out_dir, f"{out_prefix}_filtered_taxa_{rank}.tsv")
        # try:
        with open(outfile, "w") as fout:
            fout.writelines(output_filter_list)
        # except TypeError:
        #     pass
    if args.verbose > 0:
        print("==DEBUG== Combining reports: Done.")
    return True

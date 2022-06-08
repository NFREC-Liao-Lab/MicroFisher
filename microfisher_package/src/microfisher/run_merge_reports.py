import os
import sys
import warnings
import argparse

from . import taxonomy_utils
from . import merging_algorithm
from . import output_util

# report_files = ["eg_a.report.tsv", "eg_b.report.tsv"]
# report_files = ["eg1.report.tsv", "eg2.report.tsv"]

#
# def check_length_gt(length):
#     class RequiredLength(argparse.Action):
#         def __call__(self, parser, args, values, option_string=None):
#             if not length <= len(values):
#                 msg = 'argument "{f}" requires at length {length} arguments'.format(
#                     f=self.dest, length=length)
#                 raise argparse.ArgumentTypeError(msg)
#             setattr(args, self.dest, values)
#     return RequiredLength
#
#
#
# def main():
#
#     parser = argparse.ArgumentParser(description="MicroFisher: TODO XXX.", formatter_class=argparse.RawTextHelpFormatter)
#     parent_parser = argparse.ArgumentParser(description="Merge report")
#     # formatter_class=argparse.ArgumentDefaultsHelpFormatter)
#     parser.add_argument("--combine", nargs="+", required=True,
#                         help="Report file(s) to combine. minimum 2 files",
#                         metavar="report_1 report_2 [report_n ...]",
#                         action=check_length_gt(2))
#     parser.add_argument("--length", nargs="+", required=False, type=int,
#                         help="Minimum matching length used to generate reports, used in the weighting scheme",
#                         metavar="length_1 length_2 [length_n ...]",
#                         action=check_length_gt(2))
#     parser.add_argument("--mode", choices=merging_algorithm.MODE_CHOICES,
#                         default="raw",
#                         help="""Algorithm for combining results together.
#     boolean: Present or absent of the taxa.
#     raw: sum of the number of reads.
#     weighted: normalised by the total number of reads.
#     weighted_length: normalised by the total number of reads and minimum length (requires --length).
#     (probability): NOT yet implemented.
# """)
#     parser.add_argument("--filter", default=0.00001, type=float,
#                         help="filter out taxa if the proportion is less than %(default)s")
#     parser.add_argument("--out_dir", default="merged_results",
#                         help="Output folders for all results.")
#
#     args = parser.parse_args()
#     print(args)

def run(args):
    report_files = [os.path.join(args.workspace, f) for f in args.combine]
    out_dir = os.path.join(args.workspace, args.out_dir)
    threshold = args.filter
    # min_db_conut = args.min_overlap
    # length_list = args.length

    desired_ranks = ["class", "order", "family", "genus", "species"]
    parsed_data = taxonomy_utils.parse_report_files(report_files)

    try:
        os.makedirs(out_dir)
    except FileExistsError:
        pass

    if args.mode == "weighted_length":
        report_length_dict = taxonomy_utils.create_report_length(report_files, args.length)

    for rank in desired_ranks:
        data_summary = {k: taxonomy_utils.summarise_data(v, rank)
            for k, v in parsed_data.items()}

        if args.mode == "raw":
            results = merging_algorithm.a_raw(data_summary)
            percentage = merging_algorithm.calculate_percentage(results)
            output_list, output_filter_list = output_util.format_merge_two(
                results, percentage, threshold=threshold, label=["merged_reads", "percentage"])

        elif args.mode == "boolean":
            results = merging_algorithm.a_boolean(data_summary, min_db_conut=2)
            sorted_key = sorted(results, key=results.get, reverse=True)
            output_list, output_filter_list = output_util.format_merge_single(
                results, threshold=None, label="db_conut")

        elif args.mode == "weighted":
            results = merging_algorithm.a_weighted(data_summary)
            output_list, output_filter_list = output_util.format_merge_single(
                results, threshold=threshold, label="proportion")

        elif args.mode == "weighted_length":
            results = merging_algorithm.a_weighted(data_summary, report_length_dict)
            output_list, output_filter_list = output_util.format_merge_single(
                results, threshold=threshold, label="proportion")

        elif args.mode == "probability":
            sys.exit("Exit: mode not yet implemented")
            results = merging_algorithm.a_probability(data_summary)

        # output = output_header + "\n".join(output_list)
        outfile = os.path.join(out_dir, f"merged_output_taxa_{rank}.tsv")
        with open(outfile, "w") as fout:
            fout.write(output_list)

        outfile = os.path.join(out_dir, f"merged_output_filtered_taxa_{rank}.tsv")
        try:
            with open(outfile, "w") as fout:
                fout.write(output_filter_list)
        except TypeError:
            pass
    return True

import sys
import warnings
import argparse
import taxonomy_utils
import merging_algorithm

# report_files = ["eg_a.report.tsv", "eg_b.report.tsv"]
# report_files = ["eg1.report.tsv", "eg2.report.tsv"]


def check_length_gt(length):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not length <= len(values):
                msg = 'argument "{f}" requires at length {length} arguments'.format(
                    f=self.dest, length=length)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength

def format_output(results, threshold, label="proportion"):

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


def main():

    parser = argparse.ArgumentParser(description="MicroFisher: TODO XXX.", formatter_class=argparse.RawTextHelpFormatter)
    parent_parser = argparse.ArgumentParser(description="Merge report")
    # formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--combine", nargs="+", required=True,
                        help="Report file(s) to combine. minimum 2 files",
                        metavar="report_1 report_2 [report_n ...]",
                        action=check_length_gt(2))
    parser.add_argument("--length", nargs="+", required=False, type=int,
                        help="Minimum matching length used to generate reports, used in the weighting scheme",
                        metavar="length_1 length_2 [length_n ...]",
                        action=check_length_gt(2))
    parser.add_argument("--mode", choices=merging_algorithm.MODE_CHOICES,
                        default="raw",
                        help="""Algorithm for combining results together.
    boolean: Present or absent of the taxa.
    raw: sum of the number of reads.
    weighted: normalised by the total number of reads.
    weighted_length: normalised by the total number of reads and minimum length (requires --length).
    (probability): NOT yet implemented.
""")
    parser.add_argument("--filter", default=0.00001, type=float,
                        help="filter out taxa if the proportion is less than %(default)s")

    args = parser.parse_args()
    print(args)
    report_files = args.combine
    threshold = args.filter
    # min_db_conut = args.min_overlap
    # length_list = args.length
    desired_ranks = ["class", "order", "family", "genus", "species"]

    report_length_dict = taxonomy_utils.create_report_length(
        report_files, args.length)
    parsed_data = taxonomy_utils.parse_report_files(report_files)

    for rank in desired_ranks:
        data_summary = {k: taxonomy_utils.summarise_data(
            v, rank) for k, v in parsed_data.items()}

        if args.mode == "raw":
            results = merging_algorithm.a_raw(data_summary)
            percentage = merging_algorithm.calculate_percentage(results)
            sorted_key = sorted(results, key=results.get, reverse=True)

            output_header = "taxID\tname\tmerged_reads\tpercentage\n"
            output_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}"
                           for k in sorted_key]
            output_list = output_header + "\n".join(output_list)
            # filtered_key = [k for k,v in percentage.items() if v>threshold]
            filtered_key = [k for k in sorted_key if percentage[k] > threshold]
            output_filter_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}"
                                  for k in filtered_key]
            output_filter_list = output_header + "\n".join(output_filter_list)

        elif args.mode == "boolean":
            results = merging_algorithm.a_boolean(data_summary, min_db_conut=2)
            sorted_key = sorted(results, key=results.get, reverse=True)

            # output_header = "taxID\tname\tdb_counts\n"
            # output_list = [f"{k[0]}\t{k[1]}\t{results[k]}" for k in sorted_key]
            output_list, output_filter_list = format_output(results, threshold=None, label="db_conut")

        elif args.mode == "weighted":
            results = merging_algorithm.a_weighted(data_summary)
            output_list, output_filter_list = format_output(results, threshold=threshold, label="proportion")

        elif args.mode == "weighted_length":
            results = merging_algorithm.a_weighted(
                data_summary, report_length_dict)
            output_list, output_filter_list = format_output(results, threshold=threshold, label="proportion")

        elif args.mode == "probability":
            sys.exit("Exit: mode not yet implemented")
            results = merging_algorithm.a_probability(data_summary)

        # output = output_header + "\n".join(output_list)
        with open(f"merged_output_taxa_{rank}.tsv", "w") as fout:
            fout.write(output_list)

        # output_filter = output_header + "\n".join(output_filter_list)
        try:
            with open(f"merged_output_filtered_taxa_{rank}.tsv", "w") as fout:
                fout.write(output_filter_list)
        except TypeError:
            pass

if __name__ == "__main__":
    main()

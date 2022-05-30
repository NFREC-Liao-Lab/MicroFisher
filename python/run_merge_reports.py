import sys
import warnings
import argparse
import taxanomy_utils
import merging_algorithm

# report_files = ["eg1.report.tsv", "eg2.report.tsv"]

def main():

    parser = argparse.ArgumentParser(description="MicroFisher: TODO XXX.")
    parent_parser = argparse.ArgumentParser(description="Merge report")
    # formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--combine", nargs="+", required=True,
                       help="Report file(s) to combine. minimum 2 files",
                       metavar="report_1 report_2 [report_n ...]")
    parser.add_argument("--mode", choices=merging_algorithm.MODE_CHOICES,
                        default="raw",
                        help="Algorithm for combining results together. (probability) model not yet implemented.")
    parser.add_argument("--filter", default=0.00001,
                        help="filter out taxa if the proportion is less than %(default)s")

    args = parser.parse_args()
    print(args)
    report_files = args.combine
    threshold = args.filter
    desired_ranks=["class", "order", "family", "genus", "species"]

    parsed_data = taxanomy_utils.parse_report_files(report_files)
    for rank in desired_ranks:
        data_summary = {k:taxanomy_utils.summarise_data(v, rank) for k, v in parsed_data.items()}

        if args.mode == "raw":
            results = merging_algorithm.a_raw(data_summary)
            percentage = merging_algorithm.calculate_percentage(results)
            sorted_key = sorted(results, key=results.get, reverse=True)
            output_header = "taxID\tname\tmerged_reads\tpercentage\n"
            output_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}" for k in sorted_key]

            # filtered_key = [k for k,v in percentage.items() if v>threshold]
            filtered_key = [k for k in sorted_key if percentage[k]>threshold]
            output_filter_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}" for k in filtered_key]

        elif args.mode == "weighted":
            results = merging_algorithm.a_weighted(data_summary)
            sorted_key = sorted(results, key=results.get, reverse=True)
            output_header = "taxID\tname\tproportion\n"
            output_list = [f"{k[0]}\t{k[1]}\t{results[k]:.6f}" for k in sorted_key]

            filtered_key = [k for k in sorted_key if results[k]>threshold]
            output_filter_list = [f"{k[0]}\t{k[1]}\t{results[k]:.6f}" for k in filtered_key]

        else:
            # results = merging_algorithm.a_probability(data_summary)
            sys.exit("Exit: mode not yet implemented")

        output = output_header+"\n".join(output_list)
        with open(f"merged_output_taxa_{rank}.tsv", "w") as fout:
            fout.write(output)

        output_filter = output_header+"\n".join(output_filter_list)
        with open(f"merged_output_filtered_taxa_{rank}.tsv", "w") as fout:
            fout.write(output_filter)


if __name__ == "__main__":
    main()

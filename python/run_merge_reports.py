import argparse
import taxanomy_utils

# report_files = ["eg1.report.tsv", "eg2.report.tsv"]
def main():

    parser = argparse.ArgumentParser(description="MicroFisher: TODO XXX.")
    parent_parser = argparse.ArgumentParser(description="Merge report")
    # formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--combine", nargs="+", required=True,
                       help="Report file(s) to combine. minimum 2 files",
                       metavar="report_1 report_2 [report_n ...]")


    args = parser.parse_args()
    print(args)
    report_files = args.combine
    desired_ranks=["class", "order", "family", "genus", "species"]
    parsed_data = taxanomy_utils.parse_report_files(report_files)
    for rank in desired_ranks:
        data_summary = {k:taxanomy_utils.summarise_data(v, rank) for k, v in parsed_data.items()}
        data_normalised = {k:taxanomy_utils.normalise_data(v) for k, v in data_summary.items()}
    #
        results = dict()
        for data_each in data_summary.values():
            for k, v in data_each.items():
                try:
                    results[k] += v
                except KeyError as e:
                    results[k] = v
    #
    #
        # {k:v/total for k, v in results.items()}
    #
        total = sum(results.values())
        percentage = {k:v/total for k, v in results.items()}
        sorted_key = sorted(results, key=results.get, reverse=True)
        output_header = "taxID\tname\tmerged_reads\tpercentage\n"
        output_list = [f"{k[0]}\t{k[1]}\t{results[k]}\t{percentage[k]:.6f}" for k in sorted_key]
        output = output_header+"\n".join(output_list)
        with open(f"merged_output_taxa_{rank}.tsv", "w") as fout:
            fout.write(output)


if __name__ == "__main__":
    main()

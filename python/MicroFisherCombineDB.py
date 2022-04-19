import os
import argparse
import sys

from config import Config
from run_centrifuge import Centrifuge
from run_kreport import CentrifugeKReport


def main():

    parser = argparse.ArgumentParser(description="%(prog)s: combine multiple DB: TODO XXX.")
    parser.add_argument("-v", "--verbose", action='count', default=0)
    parser.add_argument("--dry", action="store_true")
    # parser.add_argument("-d", "--db", required=True,
    #                     help="Database name")
    # parser.add_argument("-m", "--min", type=int, default=120,
    #                     help="minimum matching length")
    parser.add_argument("-w", "--workspace",  # type="ascii",
                        help="path to your workspace/result files. Default '.' current location.")
    parser.add_argument("--combine", nargs="+", required=True,
                        help="Results file(s) to combine")
    parser.add_argument("--min_overlap", default=2, type=int, help="taxID present in at least (n) database. (Default: %(default)s)")
    parser.add_argument("-o", "--output", required= True, help="output file name.")
    # parser.add_argument("-1", help="forward reads")
    # parser.add_argument("-2", help="forward reads")
    # TODO(SW): Add either -1 -2 OR --prefix options later
    # parser.add_argument("--prefix", required=True,
    #                     help="used for both infiles and outfiles.\n Infiles are in the workspace [prefix_R1.fastq.gz,prefix_R2.fastq.gz].")
    # parser.add_argument("--centrifuge_path", default="",
    #                     help="path to you centrifuge program (if it is not available in $PATH)")
    # parser.add_argument("--db_path", default="",
    #                     help="path to the database (if it is not available by default). Currently append to centrifuge_path")

    args = parser.parse_args()
    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")

    min_db_conut = args.min_overlap
    infiles = [os.path.join(args.workspace, f) for f in args.combine]

    all_results = {}
    for f in infiles:
        with open(f, "r") as infile:
            id_dict = {}
            header = infile.readline()
            for line in infile.readlines():
                id = line.split("\t")[1]
                id_dict[id] = line
                # merged_result[id] = line
                try:
                    all_results[id].append(line)
                except KeyError as e:
                    all_results[id] = [line]

    combined_list = []
    for k, v in all_results.items():
        if len(v) >= min_db_conut:
            combined_list.append(k)

    with open(args.output, "w") as outfile:
        outfile.write(header)
        for k in combined_list:
            # TODO: ONLY output the first hit. Need better way to comibne them.
            outfile.write(all_results[k][0])


if __name__ == "__main__":
    main()



# python3 MicroFisherCombineDB.py -vv -w /Users/steven/NTU/Project_D1D2/D1D2_approach/stat_result_simulating_100species_80_3/ --combine simulating_100species_80_3_ITS1.short_read.reprot.tsv simulating_100species_80_3_ITS2.short_read.reprot.tsv simulating_100species_80_3_LsuD1.short_read.reprot.tsv simulating_100species_80_3_LsuD2.short_read.reprot.tsv --min_overlap 3 --output combined_result.report.tsv

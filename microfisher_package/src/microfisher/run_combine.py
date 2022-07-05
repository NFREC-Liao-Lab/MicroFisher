# Deprecated
import os

from .configuration import Config
from .run_kreport import CentrifugeKReport


def combine_reports(args):

    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")

    min_db_count = args.min_overlap
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
                except KeyError:
                    all_results[id] = [line]

    combined_list = []
    for k, v in all_results.items():
        if len(v) >= min_db_count:
            combined_list.append(k)

    with open(args.output, "w") as outfile:
        outfile.write(header)
        for k in combined_list:
            # TODO: ONLY output the first hit. Need better way to combine them.
            outfile.write(all_results[k][0])

    config = Config(args)
    config.out_report = args.output
    cent_kreport = CentrifugeKReport(config)
    if not args.dry:
        output = cent_kreport.run()

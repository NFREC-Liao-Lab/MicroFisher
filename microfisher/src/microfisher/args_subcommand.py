import os
import sys
from . import run_merge_reports
from . import run_init_db
from . import merging_algorithm
from . import taxonomy
from .configuration import Config
from .run_centrifuge import Centrifuge
from .run_kreport import CentrifugeKReport


def init_db(args):
    if not args.dry:
        is_complete = run_init_db.init_setup_db(output_dir=args.db_loc)
        if is_complete:
            print("==DEBUG== subcommand init_db completed.\n"
                  f"Prebuild database available at: {args.db_loc}\n")

def parse_output_dir(args):
    # args.workspace = os.path.realpath(args.workspace)
    try:
        if args.out_dir != os.path.realpath(args.out_dir):
            args.out_dir = os.path.join(os.path.realpath(args.workspace), args.out_dir)
        os.makedirs(args.out_dir, exist_ok=True)
    except PermissionError:
        print(f"==Error== Permission denied to create output directory: {args.out_dir}")
        sys.exit(-1)
    return args


def search_db(args):
    # print('((%s))' % args.search)
    args = parse_output_dir(args)
    config = Config(args)
    centrifuge = Centrifuge(config)
    cent_kreport = CentrifugeKReport(config)
    is_complete_k = False
    if not args.dry:
        is_complete = centrifuge.run()
        if is_complete:
            is_complete_k = cent_kreport.run()
        if is_complete_k:
            print("==DEBUG== subcommand search completed.\t"
                  f"Results at {config.out_report}.\n")
    return config.out_report


def combine_reports(args):
    if not args.dry:
        args = parse_output_dir(args)
        is_complete = run_merge_reports.run(args)
        if is_complete:
            print("==DEBUG== subcommand combine completed.\t"
                  f"Merged results at {args.out_dir} folder.\n")


def preset_algorithm(args):
    args.min = args.min
    args.filter = merging_algorithm.FILTER_DEFAULT
    args.mode = merging_algorithm.MODE_CHOICES[0]
    args.desired_ranks = taxonomy.DESIRED_RANKS
    args.cent_length = None
    args.db_length = None
    args.include_all = True
    combine_list = []
    if not args.dry:
        for db in merging_algorithm.DB_LIST[args.preset_db]:
            args.db = db
            outfile = search_db(args)
            combine_list.append(outfile)
        args.combine = combine_list
        combine_reports(args)

from . import run_merge_reports
from . import run_init_setup
from . import merging_algorithm
from .configuration import Config
from .run_centrifuge import Centrifuge
from .run_kreport import CentrifugeKReport


def init_db(args):
    if not args.dry:
        is_complete = run_init_setup.init_setup_db(output_dir=args.db_loc)
        if is_complete:
            print(f"==DEBUG== subcommand init_db completed.\tPrebuild database at {args.db_loc}.\n")


def search_db(args):
    # print('((%s))' % args.search)
    config = Config(args)
    centrifuge = Centrifuge(config)
    cent_kreport = CentrifugeKReport(config)
    if not args.dry:
        is_complete = centrifuge.run()
        if is_complete:
            is_complete_k = cent_kreport.run()
        if is_complete_k:
            print(f"==DEBUG== subcommand search completed.\tResults at {config.out_report}.\n")
    return config.out_report


def combine_reports(args):
    if not args.dry:
        is_complete = run_merge_reports.run(args)
        if is_complete:
            print(f"==DEBUG== subcommand combine completed.\tMerged results at {args.out_dir} folder.\n")


def preset_algorithm(args):
    args.min = 120
    args.filter = merging_algorithm.FILTER_DEFAULT
    args.mode = merging_algorithm.MODE_CHOICES[0]
    args.cent_length = None
    args.db_length = None
    combine_list = []
    if not args.dry:
        for db in merging_algorithm.DB_LIST[args.preset_db]:
            args.db = db
            outfile = search_db(args)
            combine_list.append(outfile)
        args.combine = combine_list
        combine_reports(args)

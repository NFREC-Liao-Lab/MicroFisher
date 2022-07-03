from . import __version__
import argparse
import sys

# from . import run_combine
from . import run_merge_reports
from . import run_init_setup
from . import merging_algorithm
from .configuration import Config
from .run_centrifuge import Centrifuge
from .run_kreport import CentrifugeKReport


def search_db(args):
    # print('((%s))' % args.search)
    config = Config(args)
    centrifuge = Centrifuge(config)
    cent_kreport = CentrifugeKReport(config)
    if not args.dry:
        is_complete = centrifuge.run()
        if is_complete:
            output = cent_kreport.run()
    if args.verbose > 0:
        print(f"\n==DEBUG== Done\tReport at {config.out_report}")
    return config.out_report


def combine_reports(args):
    is_complete = run_merge_reports.run(args)
    if args.verbose > 0 and is_complete:
        print(f"\n==DEBUG== Done\tMerged results at {args.out_dir}")


DB_LIST = {"ITS+LSU": ["ITS1_fisher", "ITS2_fisher", "LSU_D1_fisher_new", "LSU_D2_fisher_new"],
           "ITS": ["ITS1_fisher", "ITS2_fisher"],
           "LSU": ["LSU_D1_fisher_new", "LSU_D2_fisher_new"]}
FILTER_DEFAULT = 0.00001
def preset_algorithm(args):
    args.min=120
    combine_list = []
    for db in DB_LIST[args.preset_db]:
        args.db = db
        outfile = search_db(args)
        combine_list.append(outfile)
    args.combine = combine_list
    args.filter = FILTER_DEFAULT
    args.mode = "weighted_centlength_dblength"
    args.cent_length = None
    args.db_length = None
    combine_reports(args)

def init_db(args):
    is_complete = run_init_setup.init_setup_db(output_dir=args.db_loc)
    if args.verbose > 0 and is_complete:
        print(f"\n==DEBUG== Done\tPrebuild database at {args.db_loc}")


def check_length_gt(length):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not length <= len(values):
                msg = 'argument "{f}" requires at length {length} arguments'.format(
                    f=self.dest, length=length)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength



def add_args_input_group(parser):
    group_input = parser.add_argument_group("input")
    arg_infiles = group_input.add_mutually_exclusive_group(required=True)
    arg_infiles.add_argument("--prefix", 
                    help="used for both infiles and outfiles.\n Infiles are in the workspace [prefix_R1.fastq.gz, prefix_R2.fastq.gz].")
    arg_infiles.add_argument("--paired", nargs=2, 
                    help="Requires two paired-end files, f1_*.gz, f2_*.gz. Infiles are in the workspace.")
    arg_infiles.add_argument("--single", 
                    help="Single endfiles, f_single_*.gz. Infile is in the workspace.")
    return arg_infiles


def add_args_output_group(parser):
    group_output = parser.add_argument_group("output")
    group_output.add_argument("--out_dir", default="merged_results",
                    help="Output folders for all results.")
    return group_output

def add_args_centrifuge(parser, full_config=True):
    group_centrifuge = parser.add_argument_group("centrifuge")
    group_centrifuge.add_argument("--centrifuge_path", default="",
                          help="path to you centrifuge program (if it is not available in $PATH)")
    group_centrifuge.add_argument("--db_path", default="",
                          help="path to the database (if it is not available by default). Currently append to centrifuge_path")
    group_centrifuge.add_argument("--threads", type=int, default=1,
                    help="number of alignment threads to launch. (%(default)s)")

    if full_config:
        group_centrifuge.add_argument("-m", "--min", type=int, default=120,
                            help="minimum matching length (Default: %(default)s)")
    return group_centrifuge

def main():

    parser = argparse.ArgumentParser(
        prog="MicroFisher", description="%(prog)s: TODO XXX.")
    parent_parser = argparse.ArgumentParser(
        description="Parent parser.", add_help=False)
    # formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parent_parser.add_argument("-v", "--verbose", action='count', default=0,
                               help="Verbose mode (allow multiples)")
    parent_parser.add_argument("--dry", action="store_true", help="Dry run")
    parent_parser.add_argument('--version', action='version',
                               version='%(prog)s {version}'.format(version=__version__))
    parent_parser.add_argument("-w", "--workspace",  default=".",
                               help="path to your workspace/dataset. Default '.' current location.")
    # parser.set_defaults(func=help_message)

    subparsers = parser.add_subparsers(title="subcommand", dest='subcommand',
                                       description="Collection of %(prog)s functions.", help="See additional help 'subcommand --help'")

    p_init = subparsers.add_parser('init_db', parents=[parent_parser],
                                   help='Initialise centrifuge with prebuild databases.', conflict_handler='resolve')
    p_search = subparsers.add_parser('search', parents=[parent_parser],
                                     help='Search with centrifuge', conflict_handler='resolve')
    p_combine = subparsers.add_parser('combine', parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter,
                                      help='Combine results', conflict_handler='resolve')
    p_full = subparsers.add_parser('preset', parents=[p_combine], formatter_class=argparse.RawTextHelpFormatter,
                                      help='Preset pipeline', conflict_handler='resolve')

    p_init.add_argument("--db_loc", default="default_db", required=False,
                        help="Location to store the default centrifuge databases. (Default: ./%(default)s)")
    p_init.set_defaults(func=init_db)

 
    p_search.add_argument("-d", "--db", required=True, help="Database name")
    add_args_input_group(p_search)   
    add_args_centrifuge(p_search)
    p_search.set_defaults(func=search_db)

    # parser_kreport_only = subparsers.add_parser('kk' , help='kk results')
    # p_combine.set_defaults(subcommand='combine')
    p_combine.add_argument("--combine", nargs="+", required=True,
                           help="Results file(s) to combine. minimum 2 files",
                           metavar="result1 result2 [result_n ...]",
                           action=check_length_gt(2))
    add_args_output_group(p_combine)
    group_algorithm = p_combine.add_argument_group("combining")
    group_algorithm.add_argument("--mode", choices=merging_algorithm.MODE_CHOICES,
                           default="weighted_centlength_dblen",
                           help="""Algorithm for combining results together.
    - boolean: Present or absent of the taxa (optional: --min_overlap).
    - raw: sum of the number of reads.
    - weighted: normalised by the total number of reads.
    - weighted_centlength: normalised by the total number of reads and minimum 
        length used in centrifuge (--cent_length).
    - weighted_centlength_dblen: normalised by the total number of reads, and 
        minimum length used in centrifuge (--cent_length), and the average 
        length of the database (--db_length).
""")
# (probability): NOT yet implemented.
    group_algorithm.add_argument("--filter", default=FILTER_DEFAULT, type=float,
                                 help="filter out taxa if the proportion is less than %(default)s")
    group_algorithm.add_argument("--min_overlap", default=2, required=False, type=int,
                                 help="taxID present in at least (n) database, used in the 'boolean' scheme. (Default: %(default)s)")
    group_algorithm.add_argument("--cent_length", nargs="+", required=False, type=int,
                                 help="Minimum matching length used to generate reports in centrifuge, used in the 'weighted_length' scheme",
                                 metavar="length_1 length_2 [length_n ...]",
                                 action=check_length_gt(2))
    group_algorithm.add_argument("--db_length", nargs="+", required=False, type=float,
                                 help="Average length in the database, used in the 'weighted_length' scheme",
                                 metavar="length_1 length_2 [length_n ...]",
                                 action=check_length_gt(2))

    # p_combine.add_argument("--combine_db", dest="db", required=True,
    #                        help="Combined database name")
    p_combine.set_defaults(func=combine_reports)



    p_full.add_argument("--preset_db", choices=["ITS+LSU", "ITS", "LSU"], default="ITS+LSU")
    # p_full.add_argument("--prefix", required=True,
    #                       help="used for both infiles and outfiles.\n Infiles are in the workspace [prefix_R1.fastq.gz, prefix_R2.fastq.gz].")
    # p_full.add_argument("--centrifuge_path", default="",
    #                       help="path to you centrifuge program (if it is not available in $PATH)")
    # p_full.add_argument("--db_path", default="",
    #                       help="path to the database (if it is not available by default). Currently append to centrifuge_path")
    # group_input = p_full.add_argument_group("input")
    add_args_input_group(p_full)
    add_args_output_group(p_full)
    add_args_centrifuge(p_full, full_config=False)
    p_full.set_defaults(func=preset_algorithm)



    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_help()
        sys.exit(-1)

    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")
    args.func(args)


if __name__ == "__main__":
    main()

import os
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


def merge_reports(args):
    is_complete = run_merge_reports.run(args)
    if args.verbose > 0 and is_complete :
        print(f"\n==DEBUG== Done\tMerged results at {args.out_dir}")

def init_db(args):
    is_complete = run_init_setup.init_setup_db(output_dir = args.db_loc)
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

from . import __version__
def main():

    parser = argparse.ArgumentParser(prog= "MicroFisher", description="%(prog)s: TODO XXX.")
    parent_parser = argparse.ArgumentParser(
        description="Parent parser.", add_help=False)
    # formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parent_parser.add_argument("-v", "--verbose", action='count', default=0,
                               help="Verbose mode (allow multiples)")
    parent_parser.add_argument("--dry", action="store_true", help="Dry run")
    parent_parser.add_argument("-w", "--workspace",  default=".",
                               help="path to your workspace/dataset. Default '.' current location.")
    parent_parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    # parser.set_defaults(func=help_message)

    subparsers = parser.add_subparsers(title="subcommand", dest='subcommand',
                                       description="Collection of %(prog)s functions.", help="See additional help 'subcommand --help'")

    p_init = subparsers.add_parser('init_db', parents=[parent_parser],
                                   help='Initialise centrifuge with prebuild databases.', conflict_handler='resolve')
    p_search = subparsers.add_parser('search', parents=[parent_parser],
                                     help='Search with centrifuge', conflict_handler='resolve')
    p_combine = subparsers.add_parser('combine', parents=[parent_parser], formatter_class=argparse.RawTextHelpFormatter,
                                      help='Combine results', conflict_handler='resolve')

    p_init.add_argument("--db_loc", default= "default_db", required=False,
                              help="Location to store the default centrifuge databases. (Default: ./%(default)s)")
    p_init.set_defaults(func=init_db)


    # parser.add_argument("-1", help="forward reads")
    # parser.add_argument("-2", help="forward reads")
    # TODO(SW): Add either -1 -2 OR --prefix options later
    # p_search.set_defaults(subcommand='search')
    p_search.add_argument("--prefix", required=True,
                          help="used for both infiles and outfiles.\n Infiles are in the workspace [prefix_R1.fastq.gz,prefix_R2.fastq.gz].")
    p_search.add_argument("-d", "--db", required=True, help="Database name")
    p_search.add_argument("--centrifuge_path", default="",
                          help="path to you centrifuge program (if it is not available in $PATH)")
    p_search.add_argument("--db_path", default="",
                          help="path to the database (if it is not available by default). Currently append to centrifuge_path")
    p_search.add_argument("-m", "--min", type=int, default=120,
                          help="minimum matching length (Default: %(default)s)")
    p_search.set_defaults(func=search_db)


    # parser_kreport_only = subparsers.add_parser('kk', help='kk results')
    # p_combine.set_defaults(subcommand='combine')
    p_combine.add_argument("--combine", nargs="+", required=True,
                           help="Results file(s) to combine. minimum 2 files",
                           metavar="result1 result2 [result_n ...]",
                           action=check_length_gt(2))
    p_combine.add_argument("--out_dir", default="merged_results",
                           help="Output folders for all results.")
    p_combine.add_argument("--mode", choices=merging_algorithm.MODE_CHOICES,
                           default="raw",
                           help="""Algorithm for combining results together.
    boolean: Present or absent of the taxa (optional: --min_overlap).
    raw: sum of the number of reads.
    weighted: normalised by the total number of reads.
    weighted_length: normalised by the total number of reads and minimum length (requires --length).
    (probability): NOT yet implemented.
""")
    p_combine.add_argument("--filter", default=0.00001, type=float,
                           help="filter out taxa if the proportion is less than %(default)s")
    p_combine.add_argument("--min_overlap", default=2, required=False, type=int,
                           help="taxID present in at least (n) database, used in the 'boolean' scheme. (Default: %(default)s)")
    p_combine.add_argument("--length", nargs="+", required=False, type=int,
                           help="Minimum matching length used to generate reports, used in the 'weighted_length' scheme",
                           metavar="length_1 length_2 [length_n ...]",
                           action=check_length_gt(2))

    # p_combine.add_argument("--combine_db", dest="db", required=True,
    #                        help="Combined database name")
    p_combine.set_defaults(func=merge_reports)

    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_help()
        sys.exit(-1)

    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")
    args.func(args)


if __name__ == "__main__":
    main()

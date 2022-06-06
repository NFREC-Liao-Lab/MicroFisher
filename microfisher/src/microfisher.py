import os
import argparse
import sys

import run_centrifuge
from config import Config
# from run_centrifuge import Centrifuge
# from run_kreport import CentrifugeKReport
# from run_combine import combine
import run_combine
import run_init_setup


def search_db(args):
    # print('((%s))' % args.search)
    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")
    config = Config(args)
    centrifuge = Centrifuge(config)
    cent_kreport = CentrifugeKReport(config)
    if not args.dry:
        is_complete = centrifuge.run()
        if is_complete:
            output = cent_kreport.run()


def combine_results(args):
    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")
    run_combine.combine(args)

def init_db(args):
    run_init_setup.init_setup_db(output_dir = args.db_loc)


def check_length_gt(length):
    class RequiredLength(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not length <= len(values):
                msg = 'argument "{f}" requires at length {length} arguments'.format(
                    f=self.dest, length=length)
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequiredLength


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
    parent_parser.add_argument("--centrifuge_path", default="",
                               help="path to you centrifuge program (if it is not available in $PATH)")
    parent_parser.add_argument("--db_path", default="",
                               help="path to the database (if it is not available by default). Currently append to centrifuge_path")
    # parser.set_defaults(func=help_message)

    subparsers = parser.add_subparsers(title="subcommand", dest='subcommand',
                                       description="Collection of %(prog)s functions.", help="See additional help 'subcommand --help'")

    p_init = subparsers.add_parser('init_db', #parents=[parent_parser],
                                   help='Initialise centrifuge with prebuild databases.', conflict_handler='resolve')
    p_search = subparsers.add_parser('search', parents=[parent_parser],
                                     help='Search with centrifuge', conflict_handler='resolve')
    p_combine = subparsers.add_parser('combine', parents=[parent_parser],
                                      help='Combine results', conflict_handler='resolve')

    p_init.add_argument("--db_loc", default= "default_db", required=False,
                              help="Location to store the default centrifuge databases. (Default: ./%(default)s)")
    p_init.set_defaults(func=init_db)
    # parser.add_argument("-1", help="forward reads")
    # parser.add_argument("-2", help="forward reads")
    # TODO(SW): Add either -1 -2 OR --prefix options later
    # p_search.set_defaults(subcommand='search')
    p_search.add_argument("-d", "--db", required=True,
                          help="Database name")
    p_search.add_argument("--prefix", required=True,
                          help="used for both infiles and outfiles.\n Infiles are in the workspace [prefix_R1.fastq.gz,prefix_R2.fastq.gz].")
    p_search.add_argument("-m", "--min", type=int, default=120,
                          help="minimum matching length (Default: %(default)s)")

    p_search.set_defaults(func=search_db)

    # parser_kreport_only = subparsers.add_parser('kk', help='kk results')
    # p_combine.set_defaults(subcommand='combine')
    p_combine.add_argument("--combine", nargs="+", required=True,
                           help="Results file(s) to combine. minimum 2 files",
                           metavar="result1 result2 [result_n ...]",
                           action=check_length_gt(2))
    p_combine.add_argument("--min_overlap", default=2, type=int,
                           help="taxID present in at least (n) database. (Default: %(default)s)")
    p_combine.add_argument("-o", "--output", required=True,
                           help="output file name.")
    p_combine.add_argument("--combine_db", dest="db", required=True,
                           help="Combined database name")
    p_combine.set_defaults(func=combine_results)

    args = parser.parse_args()
    print(f"\n==DEBUG== Arguments: {args}")
    if args.subcommand is None:
        parser.print_help()
        sys.exit(-1)

    # if args.verbose > 0:
    #     print(f"\n==DEBUG== Arguments: {args}")
    args.func(args)


if __name__ == "__main__":
    main()

import argparse
import sys 

from . import taxonomy

from . import __version__, args_subcommand, merging_algorithm


def check_is_list():
    class RequiredList(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if not isinstance(values, list):
                try:
                    values = values.lower()
                    data = values.split(",")
                    data = [d.strip() for d in data if len(d) > 0]
                    check = [d in taxonomy.FULL_RANKS for d in data]
                    if not all(check):
                        message = f"Contain invalid taxonomy rank. Input: {values}. Parsed: {data}"
                        raise argparse.ArgumentTypeError(message)
                    setattr(args, self.dest, data)
                except AttributeError:
                    message = 'argument "{f}" requires csv format, i.e. A,B,C'.format(
                        f=self.dest)
                    raise argparse.ArgumentTypeError(message)
    return RequiredList


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
    group_input = parser.add_argument_group(
        "input (only one of the following argument is allowed)")
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
    group_output.add_argument("--out_prefix", default="reports_",
                              help="Prefix for output reports.")
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
    parent_parser.add_argument("-w", "--workspace", default=".",
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
    p_init.set_defaults(func=args_subcommand.init_db)

    p_search.add_argument("-d", "--db", required=True, help="Database name")
    add_args_input_group(p_search)
    add_args_centrifuge(p_search)
    p_search.set_defaults(func=args_subcommand.search_db)

    p_combine.add_argument("--combine", nargs="+", required=True,
                           help="Results file(s) to combine. minimum 2 files",
                           metavar="result1 result2 [result_n ...]",
                           action=check_length_gt(2))
    add_args_output_group(p_combine)
    group_algorithm = p_combine.add_argument_group("combining")
    group_algorithm.add_argument("--mode", choices=merging_algorithm.MODE_CHOICES,
                                 default=merging_algorithm.MODE_CHOICES[0],
                                 help="""Algorithm for combining results together.
    - weighted: normalised by the total number of reads, and
        minimum length used in centrifuge (--cent_length), and the average
        length of the database (--db_length).
    - boolean: Present or absent of the taxa (optional: --min_overlap).
    - raw: sum of the number of reads.
    - weighted_abundance_only: normalised by the total number of reads (testing-only).
    - weighted_centlength_only: normalised by the total number of reads and minimum
        length used in centrifuge (--cent_length) (testing-only).
""")
    # (probability): NOT yet implemented.
    group_algorithm.add_argument("--ranks", default=taxonomy.DESIRED_RANKS,
                                 action=check_is_list(),
                                 dest="desired_ranks",
                                 help="Output results for these taxonomy ranks.")
    group_algorithm.add_argument("--filter", default=merging_algorithm.FILTER_DEFAULT, type=float,
                                 help="filter out taxa if the proportion is less than %(default)s")
    group_algorithm.add_argument("--min_overlap", default=1, required=False, type=int,
                                 help="taxID present in at least (n) database, used in the 'boolean' scheme. (Default: %(default)s)")
    group_algorithm.add_argument("--cent_length", nargs="+", required=False, type=int,
                                 metavar="length_1 length_2 [length_n ...]", action=check_length_gt(2),
                                 help="Minimum matching length used to generate reports in centrifuge, used in the 'weighted_length' scheme")
    group_algorithm.add_argument("--db_length", nargs="+", required=False, type=float,
                                 metavar="length_1 length_2 [length_n ...]",
                                 action=check_length_gt(2),
                                 help="Average length in the database, used in the 'weighted_length' scheme")
    group_algorithm.add_argument("--include_all", default=False, action=argparse.BooleanOptionalAction,
                                 help="Include unfiltered results.")##
    p_combine.set_defaults(func=args_subcommand.combine_reports)

    p_full.add_argument("--preset_db", choices=merging_algorithm.DB_LIST.keys(),
                        default="ITS+LSU")
    add_args_input_group(p_full)
    add_args_output_group(p_full)
    add_args_centrifuge(p_full, full_config=False)
    p_full.set_defaults(func=args_subcommand.preset_algorithm)

    args = parser.parse_args()

    if args.subcommand is None:
        parser.print_help()
        sys.exit(-1)

    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")
    args.func(args)


if __name__ == "__main__":
    main()

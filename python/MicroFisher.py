import os
import argparse
import sys

from config import Config
from run_centrifuge import Centrifuge
from run_kreport import CentrifugeKReport


def main():

    parser = argparse.ArgumentParser(description="MicroFisher: TODO XXX.")
    parser.add_argument("-v", "--verbose", action='count', default=0)
    parser.add_argument("--dry", action="store_true")
    parser.add_argument("-d", "--db", required=True,
                        help="Database name")
    parser.add_argument("-m", "--min", type=int, default=120,
                        help="minimum matching length")
    parser.add_argument("-w", "--workspace",  # type="ascii",
                        help="path to your workspace/dataset. Default '.' current location.")
    # parser.add_argument("-1", help="forward reads")
    # parser.add_argument("-2", help="forward reads")
    # TODO(SW): Add either -1 -2 OR --prefix options later
    parser.add_argument("--prefix", required=True,
                        help="used for both infiles and outfiles.\n Infiles are in the workspace [prefix_R1.fastq.gz,prefix_R2.fastq.gz].")
    parser.add_argument("--centrifuge_path", default="",
                        help="path to you centrifuge program (if it is not available in $PATH)")
    parser.add_argument("--db_path", default="",
                        help="path to the database (if it is not available by default). Currently append to centrifuge_path")

    args = parser.parse_args()
    if args.verbose > 0:
        print(f"\n==DEBUG== Arguments: {args}")
    config = Config(args)
    centrifuge = Centrifuge(config)
    cent_kreport = CentrifugeKReport(config)
    if not args.dry:
        is_complete = centrifuge.run()
        if is_complete:
            output = cent_kreport.run()


if __name__ == "__main__":
    main()

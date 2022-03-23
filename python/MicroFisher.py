import os
import argparse
import sys


from run_centrifuge import Centrifuge

def main():

    parser = argparse.ArgumentParser(description="MicroFisher: TODO XXX.")
    parser.add_argument("-v", "--verbose", action='count', default=0)

    parser.add_argument("-d", "--db", required=True, # type="ascii",
                        choices=["ITS1", "ITS2", "D1", "D2"],
                        help="Database name")
    parser.add_argument("-w", "--workspace", #type="ascii",
                        help="path to your workspace/dataset. Default '.' current location.")
    parser.add_argument("--centrifuge", default=".",
                        help="path to you centrifuge program (if it is not available in $PATH)")
    parser.add_argument("--db_path", default=".",
                        help="path to the database (if it is not available by default)")

    args = parser.parse_args()
    print(f"Arguments: {args}")
    cc = Centrifuge(args)


if __name__ == "__main__":
    main()

# python3 MicroFisher.py -vv --db ITS1 -w ~/workspace --centrifuge ~/workspace/centrifuge --db_path example/index

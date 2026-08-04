"""This script checks if a fastq file is sorted by length"""

import sys
import os
import argparse


def is_sorted(fastq_file, reverse=False):
    """Check if a fastq file is sorted by length"""

    with open(fastq_file) as f:
        last_len = 0 if not reverse else float("inf")
        for i, line in enumerate(f, start=1):
            if i % 4 == 2:
                seq_len = len(line.strip())
                if reverse:
                    if seq_len > last_len:
                        return False
                else:
                    if seq_len < last_len:
                        return False
                last_len = seq_len

    return True


def main():
    """Main function"""
    parser = argparse.ArgumentParser(
        description="Check if a fastq file is sorted by length"
    )
    parser.add_argument("fastq_file", help="fastq file to check", default=0, nargs="?")
    parser.add_argument(
        "-r",
        "--reverse",
        action="store_true",
        help="check if the file is sorted in descending (reverse) order",
    )
    args = parser.parse_args()

    if args.fastq_file not in [0, "-"] and not os.path.exists(args.fastq_file):
        print("Error: file not found")
        sys.exit(1)

    if args.fastq_file == "-":
        args.fastq_file = 0

    if is_sorted(args.fastq_file, reverse=args.reverse):
        print("The file is sorted by length")
        sys.exit(0)
    else:
        print("The file is NOT sorted by length")
        sys.exit(1)


if __name__ == "__main__":
    main()

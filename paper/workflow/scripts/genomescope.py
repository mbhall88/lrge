import subprocess
import sys
import argparse
from pathlib import Path
import shutil


def parse_args():
    parser = argparse.ArgumentParser(
        description="Estimate genome size with KMC and genomescopre2"
    )
    parser.add_argument("input", type=str, help="Input file")
    parser.add_argument("-k", "--kmer-size", type=int, default=21, help="Kmer size")
    parser.add_argument(
        "-m",
        "--min-kmer-count",
        type=int,
        default=2,
        help="exclude k-mers occurring less than <value> times [default: %(default)s]",
    )
    parser.add_argument(
        "-M",
        "--max-count",
        type=int,
        default=255,
        help="maximal value of a counter [default: %(default)s]",
    )
    parser.add_argument(
        "-p", "--ploidy", type=int, default=1, help="Ploidy [default: %(default)s]"
    )
    parser.add_argument(
        "--tmp", type=str, help="Temporary directory", default="gscope_tmp"
    )
    parser.add_argument("-l", "--log", type=str, help="Log file", default=2)
    parser.add_argument(
        "-C", "--no-cleanup", action="store_true", help="Do not cleanup"
    )

    return parser.parse_args()


def main():
    args = parse_args()

    tmpdir = Path(args.tmp).absolute()
    tmpdir.mkdir(exist_ok=True, parents=True)
    input_fq = Path(args.input).absolute()
    out_prefix = f'{input_fq.name.split(".")[0]}.kmc'

    # Run KMC
    kmc_cmd = [
        "kmc",
        f"-k{args.kmer_size}",
        f"-ci{args.min_kmer_count}",
        f"-cs{args.max_count}",
        str(input_fq),
        out_prefix,
        str(tmpdir),
    ]
    proc = subprocess.run(kmc_cmd, capture_output=True)
    if proc.returncode != 0:
        print(proc.stderr.decode())
        sys.exit(1)

    # run kmc_tools

    hist_out = out_prefix + ".histo"

    kmc_tools_cmd = [
        "kmc_tools",
        "transform",
        out_prefix,
        "histogram",
        hist_out,
        f"-cx{args.max_count}",
    ]
    proc = subprocess.run(kmc_tools_cmd, capture_output=True)
    if proc.returncode != 0:
        print(proc.stderr.decode())
        sys.exit(1)

    # Run GenomeScope2
    with open(args.log, "w") as log:
        gscope_cmd = [
            "genomescope2",
            "-i",
            hist_out,
            "-k",
            str(args.kmer_size),
            "-p",
            str(args.ploidy),
            "-o",
            str(tmpdir),
        ]
        proc = subprocess.run(gscope_cmd, stdout=log, stderr=log)
        if proc.returncode != 0:
            print(f"Error running GenomeScope2. See {args.log}")
            sys.exit(1)

    if not args.no_cleanup:
        try:
            shutil.rmtree(tmpdir)
        except PermissionError:
            print(f"Failed to remove {tmpdir}. Please remove it manually.")


if __name__ == "__main__":
    main()

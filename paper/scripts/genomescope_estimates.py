import argparse
from loguru import logger
import subprocess
import sys
from pathlib import Path
import os


default_kmer_sizes = "15,19,21,25,29,31"
default_min_copies = "1,2,5,10"


def parse_args():
    # Manually split the arguments into two parts at '--'
    if "--" in sys.argv:
        index = sys.argv.index("--")
        args_before_dash = sys.argv[1:index]  # Arguments before '--'
        extra_args = sys.argv[index + 1 :]  # Arguments after '--'
    else:
        args_before_dash = sys.argv[1:]
        extra_args = []

    # Initialize the parser
    parser = argparse.ArgumentParser(
        description="Explore the parameter space of genomescope2."
    )
    # Add all the expected arguments
    parser.add_argument(
        "-s",
        "--script",
        type=str,
        help="The script to explore the parameter space of.",
        default="genomescope.py",
    )
    parser.add_argument(
        "-k",
        "--kmer-sizes",
        type=str,
        default=default_kmer_sizes,
        help="The kmer sizes to use [default: %(default)s]",
    )
    parser.add_argument(
        "-m",
        "--min-copies",
        type=str,
        default=default_min_copies,
        help="The minimum number of copies to use [default: %(default)s]",
    )
    parser.add_argument(
        "-o", "--outdir", type=str, help="Output directory.", default="."
    )
    parser.add_argument(
        "inputs", type=str, nargs="+", help="Input files to pass to the mash."
    )

    # Parse only the arguments before '--'
    args = parser.parse_args(args_before_dash)

    # Store the extra arguments separately
    args.extra_args = extra_args

    print(args_before_dash)
    print(extra_args)

    return args


def main():
    args = parse_args()
    logger.debug(f"Arguments: {args}")

    outdir = Path(args.outdir).absolute()
    script = Path(args.script).absolute()

    inputs = [Path(f).absolute() for f in args.inputs]

    for k in args.kmer_sizes.split(","):
        for min_copies in args.min_copies.split(","):
            for fq in inputs:
                sample = fq.name.split(".")[0].split("__")[0]
                workdir = outdir / f"k{k}" / f"m{min_copies}" / sample
                workdir.mkdir(exist_ok=True, parents=True)
                logfile = workdir / f"{sample}.log"
                jobname = f"param_explore_gscope_{sample}_k{k}_m{min_copies}"
                # create a temporary output file name for the mash sketch. We don't
                # need the output, but we don't want it clogging up the directory
                # with a bunch of files.
                tmpdir = Path(os.environ.get("TMPDIR", "/tmp"))
                tmpfile = tmpdir / f"{jobname}"

                cmd = [
                    "ssubmit",
                    "-m",
                    "60g",
                    "-t",
                    "15m",
                    str(jobname),
                    f"/usr/bin/time -v -a -o {logfile} python {script} -k {k} -m {min_copies} -M 10000 -l {logfile} --tmp /tmp {" ".join(args.extra_args)} {fq}",
                ]

                logger.info(f"Running {cmd}. See {logfile}")
                proc = subprocess.run(cmd, stderr=subprocess.PIPE, cwd=workdir)
                stderr = proc.stderr.decode("utf-8")

                if proc.returncode != 0:
                    logger.error(f"Error running {cmd} (see {logfile})")
                    logger.error(stderr)
                    sys.exit(1)

    logger.success("All done!")


if __name__ == "__main__":
    main()

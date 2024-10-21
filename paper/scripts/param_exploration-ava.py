"""This script is for exploring the parameter space of the lrge-ava.py script."""

import argparse
import tempfile
import subprocess
import time
import sys
from loguru import logger
from pathlib import Path

default_n = "500,1000,5000,10000,25000,50000"
default_strategies = "rand,long"


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
        description="Explore the parameter space of the lrge-ava.py script."
    )

    # Add all the expected arguments
    parser.add_argument(
        "--script",
        type=str,
        help="The script to explore the parameter space of.",
        default="lrge-ava.py",
    )
    parser.add_argument(
        "-n",
        "--num-reads",
        type=str,
        default=default_n,
        help="The number of reads to use for genome size estimation [default: %(default)s]",
    )
    parser.add_argument(
        "-s",
        "--strategy",
        type=str,
        default=default_strategies,
        help="Read selection strategy [default: %(default)s]",
    )
    parser.add_argument(
        "-N",
        "--reruns",
        type=int,
        default=1,
        help="Number of reruns for each parameter combination.",
    )
    parser.add_argument(
        "-t", "--threads", type=int, help="Number of threads.", default=4
    )
    parser.add_argument(
        "-o", "--outdir", type=str, help="Output directory.", default="."
    )
    parser.add_argument(
        "inputs", type=str, nargs="+", help="Input files to pass to the script."
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

    jobids = {}
    for n in map(int, args.num_reads.split(",")):
        for strategy in args.strategy.split(","):
            for fq in inputs:
                for i in range(args.reruns):
                    sample = fq.name.split(".")[0].split("__")[0]
                    workdir = outdir / strategy / f"n{n}" / sample
                    workdir.mkdir(exist_ok=True, parents=True)
                    logfile = workdir / f"{sample}_N{i}.log"
                    jobname = f"param_explore_ava_{sample}_{strategy}_n{n}_N{i}"

                    cmd = [
                        "ssubmit",
                        "-m",
                        "12g",
                        "-t",
                        "45m",
                        str(jobname),
                        f"/usr/bin/time -v -a -o {logfile} python {script} -t {args.threads} -v -n {n} -s {strategy} {" ".join(args.extra_args)} {fq} > {logfile} 2>&1",
                        "--",
                        f"-c{args.threads}",
                    ]

                    if i > 0:
                        prev_jobid = jobids[(n, strategy, fq, i - 1)]
                        cmd.append(f"--dependency=afterok:{prev_jobid}")

                    logger.info(f"Running {cmd}. See {logfile}")
                    proc = subprocess.run(cmd, stderr=subprocess.PIPE, cwd=workdir)
                    stderr = proc.stderr.decode("utf-8")

                    if proc.returncode != 0:
                        logger.error(f"Error running {cmd} (see {logfile})")
                        logger.error(stderr)
                        sys.exit(1)
                    else:
                        jobid = stderr.strip().split()[-1]
                        try:
                            jobid = int(jobid)
                        except ValueError:
                            logger.error(f"Failed to extract jobid from {stderr}")
                            sys.exit(1)
                        jobids[(n, strategy, fq, i)] = jobid

                    # we don't need to repeat the long strategy as there's no randomness
                    if strategy == "long":
                        break

    logger.success("All done!")


if __name__ == "__main__":
    main()

"""This script is for exploring the parameter space of the lrge.py script."""

import argparse
import tempfile
import subprocess
import sys
from loguru import logger
from pathlib import Path

default_longest_sizes = "50,100,500,1000,5000,10000"
default_overlap_sizes = default_longest_sizes


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
        description="Explore the parameter space of the lrge.py script."
    )

    # Add all the expected arguments
    parser.add_argument(
        "-s",
        "--script",
        type=str,
        help="The script to explore the parameter space of.",
        default="lrge.py",
    )
    parser.add_argument(
        "-L",
        "--longest",
        type=str,
        default=default_longest_sizes,
        help="The sizes of the longest read set to explore [default: %(default)s]",
    )
    parser.add_argument(
        "-O",
        "--overlap",
        type=str,
        default=default_overlap_sizes,
        help="The sizes of the overlap set to explore [default: %(default)s]",
    )
    parser.add_argument("-R", "--reruns", type=int, default=1, help="Number of reruns.")
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

    return args


def main():
    args = parse_args()
    logger.debug(f"Arguments: {args}")

    outdir = Path(args.outdir).absolute()
    script = Path(args.script).absolute()

    inputs = [Path(f).absolute() for f in args.inputs]

    jobids = {}
    for longest in map(int, args.longest.split(",")):
        for overlap in map(int, args.overlap.split(",")):
            for strat in ["rand", "long"]:
                reruns = args.reruns if strat == "rand" else 1
                for fq in inputs:
                    for i in range(reruns):
                        sample = fq.name.split(".")[0].split("__")[0]
                        workdir = outdir / strat / f"L{longest}/O{overlap}/{sample}"
                        workdir.mkdir(exist_ok=True, parents=True)
                        logfile = (
                            workdir / f"{sample}_{strat}_L{longest}_O{overlap}_N{i}.log"
                        )
                        jobname = (
                            f"param_explore_{sample}_{strat}_L{longest}_O{overlap}_N{i}"
                        )

                        strat_flag = "-r" if strat == "rand" else ""

                        cmd = [
                            "ssubmit",
                            "-m",
                            "10g",
                            "-t",
                            "30m",
                            jobname,
                            f"/usr/bin/time -v -a -o {logfile} python {script} {strat_flag} -t {args.threads} -L {longest} -O {overlap} -v {" ".join(args.extra_args)} {fq} > {logfile} 2>&1",
                            "--",
                            f"-c{args.threads}",
                        ]

                        index = (sample, strat, longest, overlap, i)
                        if i > 0:
                            prev_jobid = jobids[
                                (sample, strat, longest, overlap, i - 1)
                            ]
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
                            jobids[index] = jobid

    logger.success("All done!")


if __name__ == "__main__":
    main()

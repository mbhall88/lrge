import argparse
from loguru import logger
import subprocess
import sys
from pathlib import Path
import os


default_sketch_sizes = "1000,10000,100000,1000000"
default_min_copies = "1,2,5,8,10"


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
        description="Explore the parameter space of the mash."
    )
    # Add all the expected arguments
    parser.add_argument(
        "-s",
        "--sketch-sizes",
        type=str,
        default=default_sketch_sizes,
        help="The sketch sizes to use for the mash [default: %(default)s]",
    )
    parser.add_argument(
        "-m",
        "--min-copies",
        type=str,
        default=default_min_copies,
        help="The minimum number of copies to use for the mash [default: %(default)s]",
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

    inputs = [Path(f).absolute() for f in args.inputs]

    for sketch_size in args.sketch_sizes.split(","):
        for min_copies in args.min_copies.split(","):
            for fq in inputs:
                sample = fq.name.split(".")[0].split("__")[0]
                workdir = outdir / f"s{sketch_size}" / f"m{min_copies}" / sample
                workdir.mkdir(exist_ok=True, parents=True)
                logfile = workdir / f"{sample}.log"
                jobname = f"param_explore_mash_{sample}_s{sketch_size}_m{min_copies}"
                # create a temporary output file name for the mash sketch. We don't
                # need the output, but we don't want it clogging up the directory
                # with a bunch of files.
                tmpdir = Path(os.environ.get("TMPDIR", "/tmp"))
                tmpfile = tmpdir / f"{jobname}"

                cmd = [
                    "ssubmit",
                    "-m",
                    "12g",
                    "-t",
                    "45m",
                    str(tmpfile),
                    f'/usr/bin/time -v -a -o {logfile} mash sketch -r -m {min_copies} -s {sketch_size} -o {tmpfile} {" ".join(args.extra_args)} {fq} > {logfile} 2>&1',
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

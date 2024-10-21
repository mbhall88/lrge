import argparse
from pathlib import Path
from loguru import logger
import sys
import re
import subprocess
import shutil
from collections import defaultdict
from statistics import median


def arg_parser():
    parser = argparse.ArgumentParser(description="Long Read Genome size Estimator")
    parser.add_argument("input", help="fastq file to estimate genome size from")
    parser.add_argument(
        "-n",
        "--num-reads",
        type=int,
        default=5000,
        help="Number of reads to use for genome size estimation [default: %(default)s]",
    )
    parser.add_argument(
        "-P",
        "--platform",
        default="ont",
        choices=["ont", "pb"],
        help="Sequencing platform of the reads [default: %(default)s]",
    )
    parser.add_argument(
        "-C", "--no-cleanup", action="store_true", help="Do not cleanup temp files"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="Number of threads [default: %(default)s]",
    )
    parser.add_argument(
        "-S",
        "--presorted",
        action="store_true",
        help="Input is presorted by length. Only use with longest strategy",
    )
    parser.add_argument(
        "-s",
        "--strategy",
        default="rand",
        choices=["rand", "long"],
        help="Read selection strategy [default: %(default)s]",
    )
    parser.add_argument(
        "-o",
        "--tmpdir",
        default="./.lrge",
        help="Temporary directory for intermediate files [default: %(default)s]",
    )
    parser.add_argument(
        "-l",
        "--min-read-len",
        type=int,
        default=0,
        help="Minimum read length to consider [default: %(default)s]",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        help="Verbose logging. -v for degub -vv for trace",
        action="count",
        default=0,
    )
    return parser.parse_args()


def longest_strategy_fastq(
    fastq: Path,
    output: Path,
    num_reads: int,
    presorted: bool,
) -> tuple[int, int]:
    num_lines = num_reads * 4

    cmd = f"cat {fastq}"
    if fastq.suffix == ".gz":
        cmd = "z" + cmd

    cmd_list = [cmd]

    if presorted:
        logger.info("Using presorted input")
    else:
        logger.info("Sorting fastq by read length")
        cmd_list.extend(
            [
                "paste - - - -",
                "perl -ne '@x=split m/\\t/; unshift @x, length($x[1]); print join \"\\t\",@x;'",
                "sort -n",
                "cut -f2-",
                'tr "\\t" "\\n"',
            ]
        )

    cmd_list.append(f"head -n {num_lines}")

    # First part: calculate cumulative read length
    cmd_read_len = " | ".join(
        cmd_list
        + [
            f"tee {output}",
            "paste - - - -",
            "cut -f2",
            'tr -d "\n"',
            "wc -c",
        ]
    )

    logger.debug(f"Running {cmd_read_len} for read length calculation")

    proc_read_len = subprocess.run(
        cmd_read_len, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )

    if proc_read_len.returncode != 0:
        stderr = proc_read_len.stderr.decode()
        logger.error(
            f"Failed to extract longest reads with the following error\n{stderr}"
        )
        sys.exit(1)

    cum_read_len = int(proc_read_len.stdout.decode().strip())

    # Second part: calculate number of reads
    cmd_n_reads = " | ".join(
        [
            f"cat {output}",
            "paste - - - -",
            "wc -l",
        ]
    )

    logger.debug(f"Running {cmd_n_reads} for read count")

    proc_n_reads = subprocess.run(
        cmd_n_reads, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )

    if proc_n_reads.returncode != 0:
        stderr = proc_n_reads.stderr.decode()
        logger.error(f"Failed to count reads with the following error\n{stderr}")
        sys.exit(1)

    n_reads = int(proc_n_reads.stdout.decode().strip())

    return n_reads, cum_read_len


def random_strategy_fastq(fastq: Path, output: Path, num_reads: int) -> tuple[int, int]:
    """Randomly subsample a fastq file"""
    # First, run the command to subsample and capture read length sum
    cmd = ["rasusa",  "reads", "-n", str(num_reads), "-o", str(output), str(fastq)]

    logger.debug(f"Running {cmd}")

    proc = subprocess.run(
        cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE
    )

    stderr = proc.stderr.decode()
    if proc.returncode != 0:
        logger.error(f"Failed to subsample fastq with the following error\n{stderr}")
        sys.exit(1)

    m = re.search(r"Keeping (?P<num_reads>\d+) reads", stderr)
    if m is None:
        logger.error("Failed to extract number of reads from rasusa output\n{stderr}")
        sys.exit(1)
    n_reads = int(m.group("num_reads"))

    m = re.search(r"Kept (?P<n_bases>\d+) bases", stderr)
    if m is None:
        logger.error("Failed to extract number of bases from rasusa output\n{stderr}")
        sys.exit(1)
    cum_read_len = int(m.group("n_bases"))

    return n_reads, cum_read_len


def minimap2_overlap(fastq: Path, output: Path, preset: str, threads: int):
    cmd = f"minimap2 -t {threads} -x {preset} {fastq} {fastq} > {output}"
    logger.debug(f"Running {cmd}")
    proc = subprocess.run(
        cmd,
        shell=True,
        stderr=subprocess.PIPE,
    )

    if proc.returncode != 0:
        stderr = proc.stderr.decode()
        logger.error(f"Failed to run minimap2 with the following error\n{stderr}")
        sys.exit(1)


def per_read_estimate(
    rlen: int, cum_rlen: int, n_reads: int, n_olaps: int, olap_thresh: int = 100
) -> float:
    potential_olaps = n_reads - 1
    x = potential_olaps / n_olaps
    avg_rlen = (cum_rlen - rlen) / (n_reads - 1)

    return x * (rlen + avg_rlen - 2 * olap_thresh)


def count_overlaps(paf: Path) -> tuple[dict[str, int], dict[str, int]]:
    olaps = defaultdict(int)
    seen_pairs = set()
    rlens = dict()
    with open(paf) as fd:
        for line in fd:
            qname = line.split("\t")[0]
            tname = line.split("\t")[5]
            if qname == tname:
                continue

            pair = tuple(sorted([qname, tname]))
            if pair in seen_pairs:
                continue
            else:
                seen_pairs.add(pair)

            tlen = int(line.split("\t")[6])
            qlen = int(line.split("\t")[1])

            olaps[tname] += 1
            olaps[qname] += 1
            rlens[tname] = tlen
            rlens[qname] = qlen

    return olaps, rlens


def estimate_genome_size(paf: Path, cum_read_len: int, n_reads: int) -> int:
    olaps, rlens = count_overlaps(paf)
    estimates = []

    for name in olaps:
        n_olaps = olaps[name]
        rlen = rlens[name]
        est = per_read_estimate(rlen, cum_read_len, n_reads, n_olaps)
        logger.trace(f"Estimate for {name}: {est}")
        estimates.append(est)

    return int(median(estimates))


def read_len_filter(fastq: Path, min_len: int, output: Path):
    cmd = f"seqkit seq -m {min_len} -o {output} {fastq}"
    logger.debug(f"Running {cmd}")
    proc = subprocess.run(
        cmd,
        shell=True,
        stderr=subprocess.PIPE,
    )

    if proc.returncode != 0:
        stderr = proc.stderr.decode()
        logger.error(f"Failed to filter reads with the following error\n{stderr}")
        sys.exit(1)


def main():
    args = arg_parser()

    if not args.verbose:
        logger.remove()
        logger.add(sys.stderr, level="INFO")
    elif args.verbose == 1:
        logger.remove()
        logger.add(sys.stderr, level="DEBUG")
    else:
        logger.remove()
        logger.add(sys.stderr, level="TRACE")

    logger.debug(args)

    cleanup = not args.no_cleanup
    preset = "ava-ont" if args.platform == "ont" else "ava-pb"
    prefix = Path(args.input).name.split(".")[0]
    input_fastq = Path(args.input)
    tmpdir = Path(args.tmpdir)
    tmpdir.mkdir(exist_ok=True, parents=True)
    output_fastq = tmpdir / f"{prefix}.{args.strategy}.fq"

    if args.min_read_len > 0:
        logger.debug(f"Filtering reads with length < {args.min_read_len}")
        filter_fastq = tmpdir / f"{prefix}.{args.strategy}.filtered.fq"
        read_len_filter(input_fastq, args.min_read_len, filter_fastq)
        input_fastq = filter_fastq

    if args.strategy == "long":
        n_reads, cum_read_len = longest_strategy_fastq(
            input_fastq, output_fastq, args.num_reads, args.presorted
        )
    elif args.strategy == "rand":
        n_reads, cum_read_len = random_strategy_fastq(
            input_fastq, output_fastq, args.num_reads
        )
    else:
        raise NotImplementedError(f"Strategy {args.strategy} not implemented")

    if n_reads < args.num_reads:
        logger.warning(
            f"Requested {args.num_reads} reads but only {n_reads} reads found"
        )
    else:
        logger.debug(f"Selected {n_reads} reads")

    logger.debug(f"Total length of selected reads: {cum_read_len}")

    output_paf = tmpdir / f"{prefix}.{args.strategy}.paf"

    minimap2_overlap(output_fastq, output_paf, preset, args.threads)

    estimate = estimate_genome_size(output_paf, cum_read_len, n_reads)

    if cleanup:
        logger.debug("Cleaning up")
        shutil.rmtree(tmpdir)

    logger.success(f"Estimated genome size: {estimate:,}")
    print(estimate)


if __name__ == "__main__":
    main()

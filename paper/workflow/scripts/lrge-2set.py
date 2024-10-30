import argparse
import tempfile
import subprocess
import sys
import re
from collections import defaultdict, Counter
from pathlib import Path
from statistics import median
import shutil
import math

from loguru import logger
import pyfastx


def arg_parser():
    parser = argparse.ArgumentParser(description="Long Read Genome size Estimator")
    parser.add_argument("input", help="fastq file to estimate genome size from")
    parser.add_argument(
        "-L",
        "--num-longest",
        type=int,
        default=100,
        help="Number of longest reads to use for genome size estimation [default: %(default)s]",
    )
    parser.add_argument(
        "-O",
        "--num-overlap-reads",
        type=int,
        default=10_000,
        help="Number of reads to overlap with the longest reads [default: %(default)s]",
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
        "-S", "--presorted", action="store_true", help="Input is presorted by length"
    )
    parser.add_argument(
        "-r",
        "--random",
        action="store_true",
        help="Select the overlap read set at random. Default is the next n longest reads after the longest reads",
    )
    parser.add_argument(
        "-I",
        "--infinite",
        action="store_true",
        help="Include reads with infinite estimates in the final estimate",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    return parser.parse_args()


def sort_fastq_by_length(fastq_file: str, output: str, reverse: bool = False):
    """gzip -dc in.fq.gz |
    paste - - - - |
    perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;' |
    sort -n |
    cut -f2- |
    tr "\t" "\n" > len_sorted.fq"""
    cmds = []
    if fastq_file.endswith(".gz"):
        cmds.append(f"gzip -dc {fastq_file}")
    else:
        cmds.append(f"cat {fastq_file}")

    cmds.append("paste - - - -")
    cmds.append(
        """perl -ne '@x=split m/\t/; unshift @x, length($x[1]); print join "\t",@x;'"""
    )
    if reverse:
        cmds.append("sort -nr")
    else:
        cmds.append("sort -n")
    cmds.append("cut -f2-")
    cmds.append('tr "\t" "\n"')

    cmd = " | ".join(cmds)
    cmd += f" > {output}"

    logger.debug(f"Running {cmd}")

    proc = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

    if proc.returncode != 0:
        stderr = proc.stderr.decode()
        logger.error(f"Failed to sort fastq with the following error\n{stderr}")
        sys.exit(1)


def random_subsample(
    input_file: str, output_file: str, num_reads: int
) -> tuple[int, int]:
    """Randomly subsample a fastq file"""
    cmd = ["rasusa", "reads", "-n", str(num_reads), "-o", output_file, input_file]

    logger.debug(f"Running {cmd}")

    proc = subprocess.run(cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE)

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


def extract_overlap_reads(
    fastq_file: str,
    output: str,
    num_overlap_reads: int,
    random: bool = False,
) -> tuple[float, int]:
    if random:
        num_reads, n_bases = random_subsample(fastq_file, output, num_overlap_reads)
        fastq_file = output
        fd_out = None
    else:
        fd_out = open(output, "w")

    avg_read_len = n_bases / num_reads

    return avg_read_len, num_reads


def split_fastq(
    sorted_fastq: str,
    longest_fastq: str,
    num_longest: int,
    rest_of_fastq: str,
    random: bool,
) -> dict[str, int]:
    longest_info = dict()

    if not random:
        in_longest = True
        with open(longest_fastq, "w") as fd_longest, open(
            rest_of_fastq, "w"
        ) as fd_rest:
            for i, (header, seq, qual) in enumerate(
                pyfastx.Fastq(sorted_fastq, build_index=False, full_name=True), start=1
            ):
                seqid = header.split()[0]
                header = f"@{seqid}"
                if in_longest:
                    longest_info[seqid] = len(seq)
                    print("\n".join([header, seq, "+", qual]), file=fd_longest)
                    if i >= num_longest:
                        in_longest = False
                else:
                    print("\n".join([header, seq, "+", qual]), file=fd_rest)
    else:
        _ = random_subsample(sorted_fastq, longest_fastq, num_longest)
        for name, seq, _ in pyfastx.Fastq(longest_fastq, build_index=False):
            longest_info[name] = len(seq)

        with open(rest_of_fastq, "w") as fd_rest:
            for header, seq, qual in pyfastx.Fastq(
                sorted_fastq, build_index=False, full_name=True
            ):
                name = header.split()[0]
                if name not in longest_info:
                    print("\n".join([f"@{header}", seq, "+", qual]), file=fd_rest)

    return longest_info


def align_overlaps_to_longest(
    longest_fq: str,
    overlap_fq: str,
    output: str,
    preset: str,
    threads: int,
) -> str:
    cmd = [
        "minimap2",
        "-o",
        output,
        "--dual=yes",
        "-t",
        str(threads),
        "-x",
        preset,
        longest_fq,
        overlap_fq,
    ]

    proc = subprocess.run(cmd, stderr=subprocess.PIPE)

    if proc.returncode != 0:
        stderr = proc.stderr.decode()
        logger.error(f"Failed to align longest reads to overlap reads\n{stderr}")
        sys.exit(1)


def per_read_estimate(
    read_len: int, read_len_stat: float, n_reads: int, n_ovlaps: int, ovlap_thresh: int
) -> float:
    """Estimate genome size using the formula:
    genome_len = (potential_overlaps / n_ovlaps) * (read_len + read_len_stat - 2 * ovlap_thresh)
    where:
    read_len_stat is the mean or median read length of the overlap reads
    """
    if n_ovlaps == 0:
        return float("inf")

    potential_overlaps = n_reads
    x = potential_overlaps / n_ovlaps
    genome_len = x * (read_len + read_len_stat - 2 * ovlap_thresh)
    return genome_len


def main():
    args = arg_parser()

    tmpdir = Path(tempfile.mkdtemp())

    if not args.verbose:  # debug is the default for loguru
        logger.remove()
        logger.add(sys.stderr, level="INFO")

    logger.debug(args)

    cleanup = not args.no_cleanup
    preset = "ava-ont" if args.platform == "ont" else "ava-pb"

    prefix = Path(args.input).name.split(".")[0]
    if not args.presorted and not args.random:
        logger.info("Sorting fastq file by length...")
        sorted_fastq = str(tmpdir / f"{prefix}.sorted.fq")
        sort_fastq_by_length(args.input, sorted_fastq, reverse=True)
    else:
        sorted_fastq = args.input

    num_longest = args.num_longest
    num_overlap_reads = args.num_overlap_reads
    longest_fastq = str(tmpdir / f"{prefix}.longest.fq")
    rest_of_fastq = str(tmpdir / f"{prefix}.rest_of.fq")
    overlap_fastq = str(tmpdir / f"{prefix}.overlap.fq")

    longest_info = split_fastq(
        sorted_fastq, longest_fastq, num_longest, rest_of_fastq, args.random
    )
    if len(longest_info) < num_longest:
        logger.error(f"Only {len(longest_info)} reads in the input file")
        logger.error(
            "This could be due to too few reads in the file, or reads with duplicate sequence IDs"
        )
        sys.exit(1)

    # extract the overlap reads from the sorted fastq. They are the num_overlap_reads after the num_longest
    logger.info("Extracting overlap reads...")

    read_len_stat, num_overlap_reads = extract_overlap_reads(
        rest_of_fastq,
        overlap_fastq,
        num_overlap_reads,
        args.random,
    )

    # the python API for minimap2 does not implement the all-v-all overlap presets so
    # we have to use the command line interface and batch all of the longest reads
    # together and align them to the overlap reads. When/ if we move to Rust, we can
    # use the all-v-all overlap presets and align one read at a time to the overlaps
    # This might be useful for very large datasets and allowing us to 'stop' if our
    # estimate seems to be converging
    logger.info("Aligning longest reads to overlap reads...")

    output_paf = str(tmpdir / f"{prefix}.overlap.paf")
    align_overlaps_to_longest(
        longest_fastq,
        overlap_fastq,
        output_paf,
        preset,
        args.threads,
    )

    logger.info("Parsing PAF file for overlap info...")

    overlap_count = Counter()
    seen_pairs = set()

    with open(output_paf) as paf:
        for line in paf:
            fields = line.split("\t")
            qname = fields[0]
            tname = fields[5]
            pair = tuple(sorted([qname, tname]))
            if pair in seen_pairs:
                continue
            else:
                seen_pairs.add(pair)

            overlap_count[tname] += 1

    if num_longest > len(overlap_count):
        missing = num_longest - len(overlap_count)
        logger.warning(f"{missing} of the longest reads had no overlaps")

    logger.info("Estimating genome size...")

    estimates = []

    for tname, rlen in longest_info.items():
        n_ovlaps = overlap_count[tname]
        n_reads = num_overlap_reads
        ovlap_thresh = 100

        est = per_read_estimate(rlen, read_len_stat, n_reads, n_ovlaps, ovlap_thresh)
        logger.debug(f"Read: {tname}, Estimate: {est}")
        estimates.append(est)

    genome_size = median(estimates)
    logger.debug(
        f"Estimate when including all estimates (i.e., finite and infinite): {genome_size:,}"
    )

    # take the median of the finite estimates
    if not args.infinite:
        finite_estimates = [x for x in estimates if not math.isinf(x)]
        genome_size = median(finite_estimates)

    logger.info(f"Estimated genome size: {genome_size:,}")
    print(genome_size)

    if cleanup:
        logger.info("Cleaning up temp files...")
        shutil.rmtree(tmpdir)

    logger.success("Done!")


if __name__ == "__main__":
    main()

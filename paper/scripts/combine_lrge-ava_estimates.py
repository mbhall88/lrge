from pathlib import Path
import sys


def main():
    # Read the estimates from the files
    indir = Path(sys.argv[1])
    print(
        ",".join(
            [
                "sample",
                "strategy",
                "n_reads",
                "repeat",
                "estimate",
                "memory",
                "cpu_time",
            ]
        )
    )
    for log in indir.rglob("*.log"):
        sample = log.parts[-2]
        n_reads = log.parts[-3][1:]
        strategy = log.parts[-4]
        repeat = log.name.split("_")[-1].split(".")[0][1:]
        est = None
        memory = None
        sys_time = None
        user_time = None
        with open(log) as f:
            for line in f:
                if "SUCCESS" in line:
                    next_line = f.readline()
                    est = float(next_line.strip())
                elif "Maximum resident set size (kbytes):" in line:
                    memory = int(line.strip().split()[-1])
                elif "System time (seconds):" in line:
                    sys_time = float(line.strip().split()[-1])
                elif "User time (seconds):" in line:
                    user_time = float(line.strip().split()[-1])

        if est is None:
            raise ValueError(f"Could not find the estimate in {log}")
        if memory is None:
            raise ValueError(f"Could not find the memory in {log}")
        if sys_time is None:
            raise ValueError(f"Could not find the sys time in {log}")
        if user_time is None:
            raise ValueError(f"Could not find the user time in {log}")

        cpu_time = round(sys_time + user_time, 2)

        print(
            ",".join(
                [
                    sample,
                    strategy,
                    n_reads,
                    repeat,
                    str(est),
                    str(memory),
                    str(cpu_time),
                ]
            )
        )


if __name__ == "__main__":
    main()

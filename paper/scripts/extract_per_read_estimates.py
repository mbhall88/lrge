import sys
from pathlib import Path


def extract_ava_estimates(infile):
    estimates = []
    with open(infile, "r") as f:
        for line in f:
            if "TRACE" in line and "Estimate for" in line:
                est = float(line.strip().split()[-1])
                estimates.append(est)

    return estimates


def extract_2set_estimates(infile):
    estimates = []
    with open(infile, "r") as f:
        for line in f:
            if "DEBUG" in line and "Estimate:" in line:
                est = float(line.strip().split()[-1])
                estimates.append(est)

    return estimates


def main():
    ava_searchdir = Path(sys.argv[1])
    two_set_searchdir = Path(sys.argv[2])

    print("accession,method,estimate")

    for p in ava_searchdir.rglob("*.log"):
        estimates = extract_ava_estimates(p)
        acc = p.name.split(".")[0]
        for est in estimates:
            print(f"{acc},ava,{est}")

    for p in two_set_searchdir.rglob("*.log"):
        estimates = extract_2set_estimates(p)
        acc = p.name.split(".")[0]
        for est in estimates:
            print(f"{acc},2set,{est}")


if __name__ == "__main__":
    main()

"""This script takes a TSV file that contains a list of RefSeq assembly information. 
One of the columns is the BioSample accession for the assembly. While another is the 
sequencing platform(s) used to generate the assembly. This script will fetch the run 
accession(s) and sequencing platform for each BioSample accession. The output is the same 
as the input, but a row for each run, meaning there will likely be multiple rows for each.
However, if the assembly has no runs with long reads, it will be skipped and therefore 
not included in the output. The script will also log a warning if a BioSample is expected 
to have PacBio or Oxford Nanopore data but none is found.
"""

import argparse
import requests
import pandas as pd
import sys
import logging

asm_acc_col = "Assembly Accession"
biosample_col = "Assembly BioSample Accession"
tech_col = "Assembly Sequencing Tech"
default_fields = "run_accession,instrument_platform"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_tsv", help="Input TSV file")
    parser.add_argument(
        "-B", "--biosample-col", default=biosample_col, help="BioSample column name"
    )
    parser.add_argument(
        "-T", "--tech-col", default=tech_col, help="Sequencing tech column name"
    )
    parser.add_argument(
        "-A",
        "--asm-acc-col",
        default=asm_acc_col,
        help="Assembly accession column name",
    )
    parser.add_argument(
        "-F", "--fields", default=default_fields, help="Fields to fetch from ENA"
    )
    parser.add_argument(
        "-I", "--add-illumina", action="store_true", help="Add Illumina run accessions"
    )
    return parser.parse_args()


def fetch_run_accessions(biosample, fields=default_fields):
    """Example response body is
    [
    {"run_accession":"SRR23686740","instrument_platform":"PACBIO_SMRT"}
    ]
    """
    url = f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=json&query=sample_accession={biosample}&fields={fields}"
    headers = {"Content-type": "application/x-www-form-urlencoded"}

    response = requests.get(url, headers=headers)
    response.raise_for_status()
    return response.json()


def main():
    args = parse_args()
    df = pd.read_csv(args.input_tsv, sep="\t", index_col=args.asm_acc_col)

    # setups logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    headers = df.columns.to_list()
    headers.insert(0, args.asm_acc_col)
    headers.extend(["Run Accession", "Instrument Platform"])
    print("\t".join(headers))

    new_rows = []
    for i, (asm_acc, row) in enumerate(df.iterrows(), start=1):
        logging.info(f"Processing assembly {asm_acc} ({i}/{len(df)})")
        rows_for_asm = []
        biosample = row[args.biosample_col]
        techs = row[args.tech_col].split(";")
        should_have_pacbio = "PacBio" in techs
        has_pacbio = False
        should_have_ont = "ONT" in techs
        has_ont = False
        has_no_long_reads = True
        run_accessions = fetch_run_accessions(biosample, args.fields)
        for entry in run_accessions:
            run = entry["run_accession"]
            platform = entry["instrument_platform"]
            if "PACBIO" in platform:
                has_pacbio = True
                has_no_long_reads = False
            if "OXFORD_NANOPORE" in platform:
                has_ont = True
                has_no_long_reads = False
            if (
                "OXFORD_NANOPORE" not in platform and "PACBIO" not in platform
            ) and not args.add_illumina:
                continue

            new_row = row.copy()
            new_row["Run Accession"] = run
            new_row["Instrument Platform"] = platform
            rows_for_asm.append(new_row)
            # get row as a tuple, including the index
            fields = new_row.to_list()
            fields.insert(0, asm_acc)
            str_fields = ["" if pd.isna(f) else str(f) for f in fields]
            print("\t".join(str_fields))

        if should_have_pacbio and not has_pacbio:
            logging.warning(
                f"No PacBio data found for assembly {asm_acc} with BioSample {biosample}"
            )
        if should_have_ont and not has_ont:
            logging.warning(
                f"No Oxford Nanopore data found for assembly {asm_acc} with BioSample {biosample}"
            )
        if has_no_long_reads or not rows_for_asm:
            logging.warning(
                f"No long read data found for assembly {asm_acc} with BioSample {biosample}. Skipping..."
            )
            continue

        # new_rows.extend(rows_for_asm)
        logging.info(f"Processed assembly {asm_acc} with BioSample {biosample}")
        sys.stdout.flush()

    # new_df = pd.DataFrame(new_rows)
    # new_df.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    main()

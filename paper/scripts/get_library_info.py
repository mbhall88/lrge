"""Run this script fetches the library information for each run in the input TSV file."""

import logging
import sys
from concurrent.futures import ThreadPoolExecutor

import pandas as pd
import requests

default_fields = "run_accession,library_selection,library_source,library_strategy"
threads = 8


def fetch_library_info(run, fields=default_fields):
    """Example response body is
    [
    {"run_accession":"SRR9821893","library_source":"TRANSCRIPTOMIC","library_strategy":"RNA-Seq","library_selection":"Inverse rRNA"}
    ]
    """
    url = f"https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=json&query=run_accession={run}&fields={fields}"
    headers = {"Content-type": "application/x-www-form-urlencoded"}

    response = requests.get(url, headers=headers)
    if response.status_code != 200:
        logging.error(f"Failed to fetch library info for {run}")
    return response.json()


def main():
    input_tsv = sys.argv[1]
    df = pd.read_csv(input_tsv, sep="\t", index_col="Run Accession")

    # setups logging
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    # fetch library info, in parallel
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {executor.submit(fetch_library_info, run): run for run in df.index}
        for future in futures:
            run = futures[future]
            try:
                library_info = future.result()[0]
                selection = library_info["library_selection"]
                df.loc[run, "Library Selection"] = selection
                source = library_info["library_source"]
                df.loc[run, "Library Source"] = source
                strategy = library_info["library_strategy"]
                df.loc[run, "Library Strategy"] = strategy
                logging.info(
                    f"Fetched library info for {run}: {selection}, {source}, {strategy}"
                )
            except Exception as e:
                logging.error(f"Failed to fetch library info for {run}: {e}")

    df.to_csv(sys.stdout, sep="\t")


if __name__ == "__main__":
    main()

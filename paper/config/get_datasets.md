These are the commands used to download and organise the info about the datasets we 
test and validate on.

Download summaries for all bacterial assemblies on RefSeq

```sh
datasets summary genome taxon bacteria --as-json-lines --assembly-source RefSeq --assembly-level complete --assembly-version latest --mag exclude > bacteria.jsonl
```

this ensures that we only get summaries for 'complete' assemblies

> A Complete assembly in RefSeq means that the genome is fully assembled, with no gaps, and typically represents the entire genome. All chromosomes and extra-chromosomal elements (e.g., plasmids or organelle genomes) are completely represented.

We convert this to a TSV file and extract the information of interest

```sh
FIELDS="accession,organism-name,organism-tax-id,assminfo-bioproject,assminfo-biosample-accession,assminfo-name,organism-infraspecific-strain,assminfo-biosample-strain,assminfo-sequencing-tech,assmstats-total-number-of-chromosomes,assmstats-total-sequence-len,assmstats-genome-coverage,assminfo-biosample-project-name,checkm-completeness,checkm-completeness-percentile,checkm-contamination,assminfo-release-date"
dataformat tsv genome --inputfile bacteria.jsonl --fields $FIELDS | sed 's/"//g' > bacteria.tsv
```

We extract the assemblies where the sequencing technology is ONT or PacBio. This is a little 
hacky as the sequencing technology field is essentially free text.

```sh
rg 'SMRT|PacBio|Sequel|sequel|Oxford|ONT|Revio|ION|Pacific|pacbio|Pacbio|OXFORD|Nanopore' bacteria.tsv > bacteria_lr.tsv
head -n1 bacteria.tsv | cat - bacteria_lr.tsv > tmp.tsv
mv tmp.tsv bacteria_lr.tsv
```

See `../notebooks/notepad.ipynb` for Python code that was used to deduplicate the rows 
that have the same biosample - keeping those with the highest coverage. In addition, there 
is code to standardise the sequencing technology field.

Next, we get the run accessions for each row (BioSample), along with their sequencing platform.

A simplified version of how this is done is querying the EBI with

```
$ biosample=SAMN31564381
$ curl "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&format=tsv&query=sample_accession=${biosample}&fields=run_accession,instrument_platform"
run_accession   instrument_platform
SRR22225500     ILLUMINA
SRR22225499     OXFORD_NANOPORE
```

We do this programmatically using `../scripts/fetch_run_accessions.py`
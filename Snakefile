from pathlib import Path

from snakemake.logging import logger

SHEET = "samples.tsv"  # approved sample sheet inside a dataset folder


def find_table(folder):
    tables = [f for f in folder.glob("*.tsv") if f.name != SHEET]
    if len(tables) != 1:
        raise ValueError(f"{folder}: expected one data table (.tsv), found {len(tables)}")
    return str(tables[0])


# Every folder in Input/ is a dataset; only those with an approved sheet go further
DATASETS = {d.name: find_table(d) for d in sorted(Path("Input").iterdir()) if d.is_dir()}
APPROVED = [d for d in DATASETS if Path("Input", d, SHEET).exists()]
PENDING = [d for d in DATASETS if d not in APPROVED]

for d in PENDING:
    logger.warning(f"AWAITING APPROVAL {d}: review results/{d}/samples_draft.tsv, "
                   f"then copy it to Input/{d}/{SHEET} and run again.")

wildcard_constraints:
    dataset="[^/]+",


rule all:
    input:
        expand("results/{dataset}/samples_draft.tsv", dataset=DATASETS),
        expand("results/{dataset}/{f}.tsv", dataset=APPROVED, f=["intensities", "samples"]),


rule draft_samples:
    input:
        table=lambda w: DATASETS[w.dataset],
    output:
        "results/{dataset}/samples_draft.tsv",
    script:
        "scripts/draft_samples.R"


rule clean_table:
    input:
        table=lambda w: DATASETS[w.dataset],
        sheet=f"Input/{{dataset}}/{SHEET}",
    output:
        intensities="results/{dataset}/intensities.tsv",
        samples="results/{dataset}/samples.tsv",
    script:
        "scripts/clean_table.R"

import csv
import os
import re
from pathlib import Path

from snakemake.logging import logger

SHEET = "samples.tsv"  # approved sample sheet inside a dataset folder
CRITERIA = "significance.tsv"  # significance criteria inside a dataset folder, one DE folder per row

# Using env packages 
os.environ.update(R_LIBS=os.devnull, R_LIBS_USER=os.devnull, R_LIBS_SITE=os.devnull)


def find_table(folder):
    tables = [f for f in folder.glob("*.tsv") if f.name not in (SHEET, CRITERIA)]
    if len(tables) != 1:
        raise ValueError(f"{folder}: expected one data table (.tsv), found {len(tables)}")
    return str(tables[0])

DATASETS = {d.name: find_table(d) for d in sorted(Path("Input").iterdir()) if d.is_dir()}
APPROVED = [d for d in DATASETS if Path("Input", d, SHEET).exists()]
PENDING = [d for d in DATASETS if d not in APPROVED]

for d in PENDING:
    logger.warning(f"AWAITING APPROVAL {d}: review results/{d}/samples_draft.tsv, "
                   f"then copy it to Input/{d}/{SHEET} and run again.")


def read_criteria(dataset):
    """{name: {padj, pvalue, log2fc}} from the dataset's significance sheet; empty cells are left out."""
    path = Path("Input", dataset, CRITERIA)
    if not path.exists():
        return {}
    criteria = {}
    with open(path, newline="", encoding="utf-8-sig") as f:  # tolerates Excel's BOM and CRLF
        for row in csv.DictReader(f, delimiter="\t"):
            name = row["name"].strip()
            if not re.fullmatch(r"[\w.-]+", name) or name == "all_proteins" or name in criteria:
                raise ValueError(f"{path}: '{name}' must be a unique folder name (letters, digits, . _ -), not all_proteins")
            cutoffs = {}
            for key in ("padj", "pvalue", "log2fc"):
                value = (row.get(key) or "").strip()
                if value:
                    cutoffs[key] = float(value.replace(",", "."))  # Excel with a comma decimal separator writes 0,05
            if "padj" not in cutoffs and "pvalue" not in cutoffs:
                raise ValueError(f"{path}: '{name}' needs a padj or pvalue cutoff")
            criteria[name] = cutoffs
    return criteria


SIGNIFICANCE = {d: read_criteria(d) for d in APPROVED}

for d in APPROVED:
    if not SIGNIFICANCE[d]:
        logger.warning(f"NO SIGNIFICANCE CRITERIA {d}: add Input/{d}/{CRITERIA} (see README) "
                       f"to get significant-protein lists; only results/{d}/de/all_proteins is built.")

wildcard_constraints:
    dataset="[^/]+",
    criterion="[^/]+",


rule all:
    input:
        expand("results/{dataset}/samples_draft.tsv", dataset=DATASETS),
        expand("results/{dataset}/{f}.tsv", dataset=APPROVED, f=["intensities", "samples"]),
        expand("results/{dataset}/qc_report.html", dataset=APPROVED),
        expand("results/{dataset}/de/all_proteins", dataset=APPROVED),
        [f"results/{d}/de/{c}/{f}" for d in APPROVED for c in SIGNIFICANCE[d] for f in ["summary.tsv", "de_report.html"]],


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


rule qc_report:
    input:
        table=lambda w: DATASETS[w.dataset],
        sheet=f"Input/{{dataset}}/{SHEET}",
        intensities="results/{dataset}/intensities.tsv",
        samples="results/{dataset}/samples.tsv",
    output:
        "results/{dataset}/qc_report.html",
    script:
        "scripts/qc_report.Rmd"


rule de_fit:
    input:
        intensities="results/{dataset}/intensities.tsv",
        samples="results/{dataset}/samples.tsv",
    output:
        stats=directory("results/{dataset}/de/all_proteins"),
        model="results/{dataset}/de/model.tsv",
    script:
        "scripts/de_fit.R"


rule de_significance:
    input:
        stats="results/{dataset}/de/all_proteins",
        samples="results/{dataset}/samples.tsv",
    output:
        summary="results/{dataset}/de/{criterion}/summary.tsv",
        significant=directory("results/{dataset}/de/{criterion}/significant"),
        volcano=directory("results/{dataset}/de/{criterion}/volcano"),
    params:
        criterion=lambda w: SIGNIFICANCE[w.dataset][w.criterion],
    script:
        "scripts/de_significance.R"


rule de_report:
    input:
        stats="results/{dataset}/de/all_proteins",
        model="results/{dataset}/de/model.tsv",
        intensities="results/{dataset}/intensities.tsv",
        samples="results/{dataset}/samples.tsv",
    output:
        "results/{dataset}/de/{criterion}/de_report.html",
    params:
        criterion=lambda w: SIGNIFICANCE[w.dataset][w.criterion],
    script:
        "scripts/de_report.Rmd"

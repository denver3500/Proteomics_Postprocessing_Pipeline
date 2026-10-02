import csv
import os
import re
from pathlib import Path

from snakemake.logging import logger

SHEET = "samples.tsv"  # approved sample sheet inside a dataset folder
CRITERIA = "significance.tsv"  # significance criteria inside a dataset folder, one DE folder per row
ENRICHMENT = "enrichment.tsv"  # enrichment flags inside a dataset folder: one row per gene set collection, one column per method

# Using env packages 
os.environ.update(R_LIBS=os.devnull, R_LIBS_USER=os.devnull, R_LIBS_SITE=os.devnull)


def find_table(folder):
    tables = [f for f in folder.glob("*.tsv") if f.name not in (SHEET, CRITERIA, ENRICHMENT)]
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

# Gene set collections for enrichment.tsv: name -> MSigDB "collection:subcollection" for human (HS) and mouse (MM).
# None = human sets mapped to mouse orthologs (mouse MSigDB has no KEGG). Codes: msigdbr::msigdbr_collections()
COLLECTIONS = {
    "hallmark":     {"HS": "H",                  "MM": "MH"},
    "go_bp":        {"HS": "C5:GO:BP",           "MM": "M5:GO:BP"},
    "go_cc":        {"HS": "C5:GO:CC",           "MM": "M5:GO:CC"},
    "go_mf":        {"HS": "C5:GO:MF",           "MM": "M5:GO:MF"},
    "reactome":     {"HS": "C2:CP:REACTOME",     "MM": "M2:CP:REACTOME"},
    "wikipathways": {"HS": "C2:CP:WIKIPATHWAYS", "MM": "M2:CP:WIKIPATHWAYS"},
    "kegg":         {"HS": "C2:CP:KEGG_LEGACY",  "MM": None},
}
METHODS = ("ora", "gsea")


def read_enrichment(dataset):
    """{collection: [methods]} from the dataset's enrichment sheet; any cell but empty, no, 0 or false flags a run."""
    path = Path("Input", dataset, ENRICHMENT)
    if not path.exists():
        return {}
    runs = {}
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            row = {k.strip().lower(): (v or "").strip() for k, v in row.items() if k}
            name = row.get("collection", "")
            if name not in COLLECTIONS or name in runs:
                raise ValueError(f"{path}: '{name}' must be listed once and be one of: {', '.join(COLLECTIONS)}")
            methods = [m for m in METHODS if row.get(m, "").lower() not in ("", "no", "0", "false")]
            if methods:
                runs[name] = methods
    return runs


ENRICHMENT_RUNS = {d: read_enrichment(d) for d in APPROVED}

for d in APPROVED:
    runs = ENRICHMENT_RUNS[d]
    if not runs:
        logger.info(f"No enrichment for {d}: add Input/{d}/{ENRICHMENT} (see README) to run it.")
    elif not SIGNIFICANCE[d] and any("ora" in m for m in runs.values()):
        logger.warning(f"ORA SKIPPED {d}: ORA tests the significant proteins, so it needs Input/{d}/{CRITERIA}.")
        ENRICHMENT_RUNS[d] = {c: [m for m in methods if m != "ora"] for c, methods in runs.items() if methods != ["ora"]}

# Plot scripts behind each step's figure folders (scripts/figures/<name>.R); editing one redraws it for every dataset
QC_FIGURES = ["intensities", "contaminants", "correlation", "pca", "pca_pairs", "cv"]
DE_FIGURES = ["counts", "volcano", "top_hits", "pvalues"]
GSEA_FIGURES = ["gsea", "gsea_overview"]
ORA_FIGURES = ["ora"]

wildcard_constraints:
    dataset="[^/]+",
    criterion="[^/]+",
    collection="[^/]+",


rule all:
    input:
        expand("results/{dataset}/samples_draft.tsv", dataset=DATASETS),
        expand("results/{dataset}/{f}.tsv", dataset=APPROVED, f=["intensities", "samples"]),
        expand("results/{dataset}/qc_report.html", dataset=APPROVED),
        expand("results/{dataset}/de/all_proteins", dataset=APPROVED),
        [f"results/{d}/de/{c}/{f}" for d in APPROVED for c in SIGNIFICANCE[d] for f in ["summary.tsv", "de_report.html"]],
        [f"results/{d}/enrichment/{c}/{c}_report.html" for d in APPROVED for c in ENRICHMENT_RUNS[d]],


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


rule qc:
    input:
        table=lambda w: DATASETS[w.dataset],
        sheet=f"Input/{{dataset}}/{SHEET}",
        intensities="results/{dataset}/intensities.tsv",
        samples="results/{dataset}/samples.tsv",
        figures=expand("scripts/figures/{name}.R", name=QC_FIGURES),
    output:
        directory("results/{dataset}/qc"),
    script:
        "scripts/qc.R"


rule qc_report:
    input:
        qc="results/{dataset}/qc",
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
        intensities="results/{dataset}/intensities.tsv",
        figures=expand("scripts/figures/{name}.R", name=DE_FIGURES),
    output:
        summary="results/{dataset}/de/{criterion}/summary.tsv",
        significant=directory("results/{dataset}/de/{criterion}/significant"),
        figures=directory("results/{dataset}/de/{criterion}/figures"),
    params:
        criterion=lambda w: SIGNIFICANCE[w.dataset][w.criterion],
    script:
        "scripts/de_significance.R"


rule de_report:
    input:
        summary="results/{dataset}/de/{criterion}/summary.tsv",
        significant="results/{dataset}/de/{criterion}/significant",
        figures="results/{dataset}/de/{criterion}/figures",
        model="results/{dataset}/de/model.tsv",
    output:
        "results/{dataset}/de/{criterion}/de_report.html",
    params:
        criterion=lambda w: SIGNIFICANCE[w.dataset][w.criterion],
    script:
        "scripts/de_report.Rmd"


rule gene_sets:
    input:
        intensities="results/{dataset}/intensities.tsv",
    output:
        "results/{dataset}/enrichment/{collection}/gene_sets.tsv",
    params:
        msigdb=lambda w: COLLECTIONS[w.collection],
    script:
        "scripts/gene_sets.R"


rule gsea:
    input:
        stats="results/{dataset}/de/all_proteins",
        samples="results/{dataset}/samples.tsv",
        gene_sets="results/{dataset}/enrichment/{collection}/gene_sets.tsv",
        figures=expand("scripts/figures/{name}.R", name=GSEA_FIGURES),
    output:
        directory("results/{dataset}/enrichment/{collection}/gsea"),
    script:
        "scripts/gsea.R"


rule ora:
    input:
        stats="results/{dataset}/de/all_proteins",
        samples="results/{dataset}/samples.tsv",
        gene_sets="results/{dataset}/enrichment/{collection}/gene_sets.tsv",
        figures=expand("scripts/figures/{name}.R", name=ORA_FIGURES),
    output:
        directory("results/{dataset}/enrichment/{collection}/ora/{criterion}"),
    params:
        criterion=lambda w: SIGNIFICANCE[w.dataset][w.criterion],
    script:
        "scripts/ora.R"


def enrichment_runs(w, method):
    return method in ENRICHMENT_RUNS[w.dataset][w.collection]


rule enrichment_report:
    input:
        gene_sets="results/{dataset}/enrichment/{collection}/gene_sets.tsv",
        samples="results/{dataset}/samples.tsv",
        intensities="results/{dataset}/intensities.tsv",
        gsea=lambda w: [f"results/{w.dataset}/enrichment/{w.collection}/gsea"] if enrichment_runs(w, "gsea") else [],
        ora=lambda w: [f"results/{w.dataset}/enrichment/{w.collection}/ora/{c}"
                       for c in SIGNIFICANCE[w.dataset]] if enrichment_runs(w, "ora") else [],
    output:
        "results/{dataset}/enrichment/{collection}/{collection}_report.html",
    params:
        msigdb=lambda w: COLLECTIONS[w.collection],
        criteria=lambda w: SIGNIFICANCE[w.dataset],
    script:
        "scripts/enrichment_report.Rmd"

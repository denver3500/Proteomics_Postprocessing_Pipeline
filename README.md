I. Running the pipeline 

0. Create the environment "conda env create -f environment.yaml"
   
1. Create a folder and put counts table in it:  Input/New_Dataset
2. Run "conda run -n proteomics snakemake -n" for dry-run check

3. Run "conda run -n proteomics snakemake -c4" for actual run. -c4 flag is the number of CPU cores used, adjust as nessecary

Old datasets are unchanged and are not re-computed unless .tsv files in input is changed (Samples or Counts)

II. Step 1 - Cleaning data

We take counts data (proteinLevelData.tsv and pivot_pack.tsv) and try to deduce conditions, samples, and amount of replicates. We write intensities.tsv and samples_draft.tsv that we manually check for any problems. After problems fixed we copy samples_draft.tsv, rename it to samples.tsv and put it to input (e.g. Input/Dataset/samples.tsv)

III. Step 2 - QC

For every approved dataset we write one folder, results/Dataset/qc/, with the report (qc_report.html), its tables (overview.tsv, sample_metrics.tsv, condition_metrics.tsv) and figures/.

The label column of samples.tsv is the readable name of a condition, shown in every report and plot (e.g. "M0 TRAP1 KO, DMEM" instead of M0_KO_DMEM). Use one label per condition; if empty, the condition with spaces instead of underscores is used. Files and folders always use the condition itself.

To drop a sample, write a reason in the exclude column of Input/Dataset/samples.tsv and run the pipeline again; both steps are recomputed and the report lists the excluded samples.

IV. Step 3 - Differential expression

Every pair of conditions is compared with limma (one model for the whole dataset). A contrast is named A_vs_B, where A is the condition listed first in samples.tsv; positive log2fc means higher in A. To flip a comparison, reorder the rows of Input/Dataset/samples.tsv.

What counts as significant is set per dataset in Input/Dataset/significance.tsv (tab-separated, editable in Excel). Each row becomes its own folder:

```
name                 padj   pvalue  log2fc
padj0.05_log2fc1     0.05           1
padj0.05_log2fc0.58  0.05           0.585
p0.01_log2fc1               0.01    1
```

A protein is significant when it is below every p-value cutoff given and |log2fc| >= log2fc (1 = 2-fold, 0.585 = 1.5-fold). Empty cells do not filter, but every row needs padj or pvalue. Editing a row recomputes only that folder; a deleted row leaves its old folder behind, so delete it by hand. Without significance.tsv only the full statistics are computed.

Paired samples (e.g. several tissues from the same mouse): add a block column to Input/Dataset/samples.tsv with the subject ID. Samples with the same block are treated as related; leave it empty for a sample that has no partner.

Output, results/Dataset/de/:
- model.tsv: how the model was fitted (pairing, within-block correlation)
- all_proteins/A_vs_B.tsv: statistics for every protein
- Criterion/summary.tsv: up and down counts per contrast
- Criterion/significant/A_vs_B.tsv: significant proteins only
- Criterion/figures/: counts, volcano/A_vs_B (one folder per contrast), top-hit heatmap, p-value histograms (see V. Figures)
- Criterion/de_report.html: counts and all figures in one page

V. Step 4 - Enrichment

Which gene sets to test, and how, is flagged per dataset in Input/Dataset/enrichment.tsv (tab-separated, editable in Excel). One row per collection; put an x under each method to run, leave it empty (or write no) to skip:

```
collection    ora   gsea
hallmark      x     x
go_bp         x     x
reactome            x
kegg          x
```

Collections: hallmark, go_bp, go_cc, go_mf, reactome, wikipathways, kegg (MSigDB, through the msigdbr package). The species comes from the protein names (_HUMAN, _MOUSE). Mouse uses the mouse MSigDB, except kegg, which the mouse MSigDB lacks: there human KEGG sets are mapped to mouse orthologs. To offer another collection, add it to COLLECTIONS in the Snakefile.

- GSEA: all proteins ranked by the limma t of each contrast; positive NES = the set is higher in A. Needs no significance criteria.
- ORA: are the up (or down) significant proteins over-represented in a set, compared with all detected proteins in the collection? Runs once for every row of significance.tsv.

Mixed sets: a set whose proteins go both ways (some up, some down) cancels out in GSEA and is diluted in ORA's up and down lists. So both methods also test the two directions together: GSEA on |t|, ORA on all significant proteins (direction "any"). A set is called mixed when it is enriched only in that combined test, not up or down alone, and its proteins include both directions. Figures show mixed sets in their own panel, with each bar split by how many proteins go each way.

A protein group counts as its first gene. Sets are tested when 10-500 of their genes were detected, and are enriched at padj < 0.05 (MIN_SET_SIZE, MAX_SET_SIZE, ENRICH_PADJ in scripts/common.R). The first run downloads each collection, so it needs internet. On an environment created before this step, run "conda env update -n proteomics -f environment.yaml" first.

Output, results/Dataset/enrichment/:
- enrichment_report.html: every collection in one page (counts per contrast, overview and figures)

and one folder per collection, Collection/:
- gene_sets.tsv: the gene sets used (genes detected in this dataset only)
- gsea/A_vs_B.tsv: every set with NES, p-value, padj and leading-edge genes; any_* columns for the |t| run; n_up/n_down = leading-edge proteins going each way; call = up, down or mixed. gsea/summary.tsv: enriched sets per contrast
- ora/Criterion/A_vs_B.tsv: sets sharing a significant protein, tested with the up, down and any (both together) lists, with fold enrichment, p-value, padj, those proteins' genes and n_up/n_down
- gsea/figures/, ora/Criterion/figures/: one folder per contrast, plus gsea/figures/overview (all contrasts, see VI. Figures)

VI. Figures

Every picture in the reports has its own folder with its data and script, so it can be restyled without rerunning the pipeline:

```
results/Dataset/qc/figures/pca/                          data.tsv  plot.R  pca.svg  pca.png
results/Dataset/de/Criterion/figures/volcano/A_vs_B/     data.tsv  contrast.tsv  plot.R  volcano_A_vs_B.svg  volcano_A_vs_B.png
```

QC figures: intensities, contaminants, correlation, pca, pca_pairs, cv. DE figures: counts, volcano/A_vs_B (one folder per contrast), top_hits, pvalues. Enrichment figures: gsea/figures/A_vs_B and ora/Criterion/figures/A_vs_B (one folder per contrast, both drawn by scripts/figures/enrichment.R; only enriched sets are shown, and a note says so when there are none), gsea/figures/overview.

- data.tsv (+ samples.tsv / contrast.tsv for some plots): exactly what is plotted, opens in Excel
- plot.R: a short ggplot/pheatmap script with settings at the top (size, font size, colours, labels; for volcanos, genes to always label in HIGHLIGHT)
- .svg for papers and editing in Inkscape/Illustrator, .png for slides

To make a version for a talk or a paper, copy the folder anywhere (results/ is overwritten on reruns), edit the settings in plot.R and run it inside the folder:

  conda run -n proteomics Rscript --vanilla plot.R    (--vanilla keeps packages from a personal R library out)

To change the default look for every dataset, edit the master script in scripts/figures/ (e.g. scripts/figures/volcano.R) and run the pipeline again: only the figures and reports that use that script are redrawn.

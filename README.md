I. Running the pipeline 

0. Create the environment "conda env create -f environment.yaml"
   
1. Create a folder and put counts table in it:  Input/New_Dataset
2. Run "conda run -n proteomics snakemake -n" for dry-run check

3. Run "conda run -n proteomics snakemake -c4" for actual run. -c4 flag is the number of CPU cores used, adjust as nessecary

Old datasets are unchanged and are not re-computed unless .tsv files in input is changed (Samples or Counts)

II. Step 1 - Cleaning data

We take counts data (proteinLevelData.tsv and pivot_pack.tsv) and try to deduce conditions, samples, and amount of replicates. We write intensities.tsv and samples_draft.tsv that we manually check for any problems. After problems fixed we copy samples_draft.tsv, rename it to samples.tsv and put it to input (e.g. Input/Dataset/samples.tsv)

III. Step 2 - QC

For every approved dataset we write one report, results/Dataset/qc_report.html, and its tables and figures in results/Dataset/qc/ (overview.tsv, sample_metrics.tsv, condition_metrics.tsv, figures/).

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

V. Figures

Every picture in the reports has its own folder with its data and script, so it can be restyled without rerunning the pipeline:

```
results/Dataset/qc/figures/pca/                          data.tsv  plot.R  pca.svg  pca.png
results/Dataset/de/Criterion/figures/volcano/A_vs_B/     data.tsv  contrast.tsv  plot.R  volcano_A_vs_B.svg  volcano_A_vs_B.png
```

QC figures: intensities, contaminants, correlation, pca, pca_pairs, cv. DE figures: counts, volcano/A_vs_B (one folder per contrast), top_hits, pvalues.

- data.tsv (+ samples.tsv / contrast.tsv for some plots): exactly what is plotted, opens in Excel
- plot.R: a short ggplot/pheatmap script with settings at the top (size, font size, colours, labels; for volcanos, genes to always label in HIGHLIGHT)
- .svg for papers and editing in Inkscape/Illustrator, .png for slides

To make a version for a talk or a paper, copy the folder anywhere (results/ is overwritten on reruns), edit the settings in plot.R and run it inside the folder:

  conda run -n proteomics Rscript --vanilla plot.R    (--vanilla keeps packages from a personal R library out)

To change the default look for every dataset, edit the master script in scripts/figures/ (e.g. scripts/figures/volcano.R) and run the pipeline again: only the figures and reports that use that script are redrawn.

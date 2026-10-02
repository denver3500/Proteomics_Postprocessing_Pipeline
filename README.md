I. Running the pipeline 

0. Create the environment "conda env create -f environment.yaml"
   
1. Create a folder and put counts table in it:  Input/New_Dataset
2. Run "conda run -n proteomics snakemake -n" for dry-run check

3. Run "conda run -n proteomics snakemake -c4" for actual run. -c4 flag is the number of CPU cores used, adjust as nessecary

Old datasets are unchanged and are not re-computed unless .tsv files in input is changed (Samples or Counts)

II. Step 1 - Cleaning data

We take counts data (proteinLevelData.tsv and pivot_pack.tsv) and try to deduce conditions, samples, and amount of replicates. We write intensities.tsv and samples_draft.tsv that we manually check for any problems. After problems fixed we copy samples_draft.tsv, rename it to samples.tsv and put it to input (e.g. Input/Dataset/samples.tsv)

III. Step 2 - QC

For every approved dataset we write one report: results/Dataset/qc_report.html 

To drop a sample, write a reason in the exclude column of Input/Dataset/samples.tsv and run the pipeline again; both steps are recomputed and the report lists the excluded samples.



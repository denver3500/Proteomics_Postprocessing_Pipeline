I. Running the pipeline 

1. Create a folder and put counts table in it:  Input/New_Dataset
2. Run "snakemake -n" for dry-run check
3. Run "snakemake -c4" for actual run. -c4 flag is the number of CPU cores used, adjust as nessecary

Old datasets are unchanged and are not re-computed unless .tsv files in input is changed (Samples or Counts)

II. Step 1 - Cleaning data

We take counts data (proteinLevelData.tsv and pivot_pack.tsv) and try to deduce conditions, samples, and amount of replicates. We write intensities.tsv and samples_draft.tsv that we manually check for any problems. After problems fixed we copy samples_draft.tsv, rename it to samples.tsv and put it to input (e.g. Input/Dataset/samples.tsv)



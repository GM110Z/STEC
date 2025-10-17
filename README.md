# STEC
Antiphage defence vs prophages encoding toxins 

1. **ncbi-datasets-batch.sh** : loop for ncbi datasers program to run on batch genomes
2. **padloc-parallel-run.sh**:Runs PADLOC prediction of anti-phage systems. Script found at : https://github.com/GM110Z/Phage-defence-scripts/blob/main/PADLOC-batch-runs/padloc-parallel-run.sh
3. **merge_pfam_padloc_defensefinder.py**. A variation of SantasHelper.py where only pfam (from edison.py), padloc and defensefinder are combined.
   Runs as ./merge_pfam_padloc_defensefinder.py --pfam pfam.tsv --padloc padloc.tsv  --antidefense defensefinder.tsv -o Final-R-file.tsv   (Optional argument is --amrfinder)
   Useful resources: https://github.com/GM110Z/PAPI-islands-analysis/blob/main/edison.py, https://github.com/GM110Z/PAPI-islands-analysis/tree/main/ExampleInputs/SantasHelper.py

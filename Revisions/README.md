**meta-fetch.py**: Uses GCA_xxx assemblies and entrez to retrieve metadata on geography, source, and host 
**source-quantification.py**: Uses output from meta-fetch.py to extract Biosample data and quantify the isolation source for a series of isolates

python source-quantification.py metadata_file.txt

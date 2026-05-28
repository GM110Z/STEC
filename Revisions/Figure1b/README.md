**meta-fetch.py**: Uses GCA_xxx assemblies and entrez to retrieve metadata on geography, source, and host 

run as: python meta-fetch.py list-ids.txt --email youremail@gmail.com

**source-quantification-toxin.py/source-quantification-phylogroup.py**: Uses output from meta-fetch.py to extract Biosample data and quantify the isolation source for a series of isolates

run as: python source-quantification.py metadata_file.txt

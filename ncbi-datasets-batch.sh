#!/bin/bash

input=$1
batch_size=10
counter=0

while read -r accession; do
    echo "Downloading $accession..."

    # Include GFF and protein FASTA (FAA) in the download
    datasets download genome accession "$accession" \
        --include gff3,protein \                       #edit this accordingly to what file type you need (e.g. fna, gb)
        --filename "${accession}.zip"

    unzip -o "${accession}.zip" -d "${accession}" && rm "${accession}.zip"

    echo "Done with $accession."

    ((counter++))

    if (( counter % batch_size == 0 )); then
        echo "Batch of $batch_size completed. Pausing for 1 minute to avoid overloading NCBI..."
        sleep 60
    else
        sleep 5  #  pause between requests
    fi
done < "$input"

echo "All downloads complete."

sleep 10

echo "Moving all downloaded files in the same folder named all_faa_gff'"
# Create destination folder
mkdir -p all_faa_gff

# Find and process all faa and gff3 files
find . -type f \( -name "*.faa" -o -name "*.gff" \) | while read -r file; do
    # Extract the GCF_xxx folder name
    gcf_id=$(basename "$(dirname "$file")")
    
    # Get file extension (.faa or .gff3)
    ext="${file##*.}"
    
    # Define new filename
    new_name="${gcf_id}.${ext}"
    
    # Copy and rename to central folder
    cp "$file" "all_faa_gff/${new_name}"
    
    echo "Moved: $file -> all_faa_gff/${new_name}"
done

echo "All .faa and .gff3 files have been renamed and moved to 'all_faa_gff/'"


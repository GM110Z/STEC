#!/bin/bash

# Usage: ./ncbidatasets.sh accession_list.txt
input="$1"
batch_size=10
counter=0

# File types to include (edit this line as needed)
INCLUDE_FILES="gff3,protein"  # options: genome, cds, rna, protein, gff3, gbff, fna

if [[ ! -f "$input" ]]; then
    echo "Error: Input file '$input' not found."
    exit 1
fi

while read -r accession; do
    [[ -z "$accession" ]] && continue  # skip empty lines

    echo "Downloading $accession..."

    datasets download genome accession "$accession" \
        --include "$INCLUDE_FILES" \
        --filename "${accession}.zip"

    if [[ -f "${accession}.zip" ]]; then
        unzip -o "${accession}.zip" -d "${accession}"
        rm "${accession}.zip"
        echo "Done with $accession."
    else
        echo "Warning: Failed to download or unzip $accession"
    fi

    ((counter++))

    if (( counter % batch_size == 0 )); then
        echo "Batch of $batch_size completed. Pausing for 1 minute to avoid overloading NCBI..."
        sleep 60
    else
        sleep 5
    fi
done < "$input"

echo "All downloads complete."

# Organise FAA and GFF files
echo "Moving all .faa and .gff3 files to 'all_faa_gff/'..."

mkdir -p all_faa_gff

find . -type f \( -name "*.faa" -o -name "*.gff" \) | while read -r file; do
    gcf_id=$(basename "$(dirname "$file")")
    ext="${file##*.}"
    new_name="${gcf_id}.${ext}"
    cp "$file" "all_faa_gff/${new_name}"
    echo "Moved: $file -> all_faa_gff/${new_name}"
done

echo "All .faa and .gff3 files have been renamed and moved to 'all_faa_gff/'"


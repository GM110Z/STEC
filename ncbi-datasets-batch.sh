#!/bin/bash

# Usage: ./ncbidatasets.sh accession_list.txt
input="$1"
batch_size=10
counter=0

# File types to include (edit this line as needed)
# genome -> genomic .fna; gbff -> GenBank .gbff; protein -> .faa; gff3 -> .gff/.gff3
INCLUDE_FILES="gff3,protein,genome,gbff"

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

# Organise FAA, GFF, GBFF, and FNA files
echo "Collecting .faa, .gff/.gff3, .gbff, and .fna files into 'all_seq_renamed/'..."

mkdir -p all_seq_renamed

# Find target files (handles both .gff and .gff3)
find . -type f \( -name "*.faa" -o -name "*.gff" -o -name "*.gff3" -o -name "*.gbff" -o -name "*.fna" \) | while read -r file; do
    # Try to derive a sensible ID for the file's parent accession.
    # Prefer a GCF_/GCA_ pattern if present anywhere in the path; fallback to immediate parent dir.
    path_id=$(echo "$file" | grep -oE 'G[CA]F_[0-9]+\.[0-9]+' | head -n1)
    if [[ -z "$path_id" ]]; then
        path_id=$(basename "$(dirname "$file")")
    fi

    ext="${file##*.}"

    # Normalise GFF extension to gff3 if needed
    if [[ "$ext" == "gff" ]]; then
        ext="gff3"
    fi

    new_name="${path_id}.${ext}"
    cp "$file" "all_seq_renamed/${new_name}"
    echo "Moved: $file -> all_seq_renamed/${new_name}"
done

echo "All files have been renamed and moved to 'all_seq_renamed/'."



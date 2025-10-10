#!/bin/bash
# Usage: ./ncbidatasets.sh accession_list.txt
set -euo pipefail
input="$1"
batch_size=10
counter=0
INCLUDE_FILES="gff3,protein,genome,gbff"

[[ -f "$input" ]] || { echo "Error: Input file '$input' not found." >&2; exit 1; }

# Clean a line: remove BOM/ZWNBSP/CR, trim, squeeze spaces
clean_line() {
  # remove UTF-8 BOM and zero-width no-break space, strip CR
  sed $'s/\uFEFF//g;s/\u200B//g' <<<"$1" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}

is_valid_accession() {
  [[ "$1" =~ ^G[AC]F_[0-9]+\.[0-9]+$ ]]
}

> skipped.txt

while IFS= read -r raw; do
  acc=$(clean_line "$raw")
  [[ -z "${acc}" ]] && continue
  if ! is_valid_accession "$acc"; then
    echo "Skipping invalid entry: '$raw'" | tee -a skipped.txt
    continue
  fi

  echo "Downloading ${acc}..."
  if datasets download genome accession "$acc" \
        --include "$INCLUDE_FILES" \
        --filename "${acc}.zip"; then
    if [[ -f "${acc}.zip" ]]; then
      unzip -o "${acc}.zip" -d "${acc}" >/dev/null
      rm -f "${acc}.zip"
      echo "Done with ${acc}."
    else
      echo "Warning: No ZIP produced for ${acc}" | tee -a skipped.txt
    fi
  else
    echo "Warning: datasets failed for ${acc}" | tee -a skipped.txt
  fi

  ((counter++))
  if (( counter % batch_size == 0 )); then
    echo "Batch of $batch_size completed. Pausing 60s..."
    sleep 60
  else
    sleep 5
  fi
done < "$input"

echo "All downloads complete."

echo "Collecting .faa, .gff/.gff3, .gbff, and .fna into 'all_seq_renamed/'..."
mkdir -p all_seq_renamed

# find and normalize names
while IFS= read -r -d '' file; do
  # derive accession id from path if present
  path_id=$(grep -oE 'G[CA]F_[0-9]+\.[0-9]+' <<<"$file" | head -n1)
  [[ -z "$path_id" ]] && path_id=$(basename "$(dirname "$file")")

  ext="${file##*.}"
  [[ "$ext" == "gff" ]] && ext="gff3"

  cp -f "$file" "all_seq_renamed/${path_id}.${ext}"
  echo "Moved: $file -> all_seq_renamed/${path_id}.${ext}"
done < <(find . -type f \( -name "*.faa" -o -name "*.gff" -o -name "*.gff3" -o -name "*.gbff" -o -name "*.fna" \) -print0)

echo "Done. See 'skipped.txt' for any rejects."



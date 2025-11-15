#!/bin/bash
# Usage: ./ncbidatasets.sh accession_list.txt
set -uo pipefail  # Removed set -e to avoid premature exit
input="$1"
batch_size=10
counter=0
INCLUDE_FILES="gff3,protein,genome,gbff"
[[ -f "$input" ]] || { echo "Error: Input file '$input' not found." >&2; exit 1; }
clean_line() {
  sed $'s/\uFEFF//g;s/\u200B//g' <<<"$1" | tr -d '\r' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}
is_valid_accession() {
  [[ "$1" =~ ^G[AC]F_[0-9]+\.[0-9]+$ ]]
}
> skipped.txt
while IFS= read -r raw; do
  echo "Processing raw line: '$raw'"
  acc=$(clean_line "$raw")
  echo "Cleaned accession: '$acc'"
  [[ -z "${acc}" ]] && { echo "Skipping empty line" | tee -a skipped.txt; continue; }
  if ! is_valid_accession "$acc"; then
    echo "Skipping invalid entry: '$raw'" | tee -a skipped.txt
    continue
  fi
  echo "Downloading ${acc}..."
  if datasets download genome accession "$acc" \
        --include "$INCLUDE_FILES" \
        --filename "${acc}.zip" 2> "${acc}_error.log"; then
    if [[ -f "${acc}.zip" ]]; then
      echo "ZIP file created: ${acc}.zip"
      if unzip -o "${acc}.zip" -d "${acc}" > "${acc}_unzip.log" 2>&1; then
        rm -f "${acc}.zip"
        echo "Done with ${acc}."
      else
        echo "Warning: unzip failed for ${acc}. See ${acc}_unzip.log" | tee -a skipped.txt
        cat "${acc}_unzip.log" | tee -a skipped.txt
        continue
      fi
    else
      echo "Warning: No ZIP produced for ${acc}" | tee -a skipped.txt
      cat "${acc}_error.log" | tee -a skipped.txt
      continue
    fi
  else
    echo "Warning: datasets failed for ${acc}. See ${acc}_error.log" | tee -a skipped.txt
    cat "${acc}_error.log" | tee -a skipped.txt
    continue
  fi
  ((counter++))
  echo "Completed $counter entries"
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

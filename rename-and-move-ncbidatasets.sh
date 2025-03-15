
#!/bin/bash

# Set the base directory
BASE_DIR=$1
DEST_DIR=$2

# Loop through each subfolder
for folder in "$BASE_DIR"/*; do
    if [ -d "$folder" ]; then
        basename=$(basename "$folder")  # Get folder name

        # Rename and move .faa files
        for file in "$folder"/*.faa; do
            if [ -f "$file" ]; then
                mv "$file" "$DEST_DIR/${basename}.faa"
            fi
        done

        # Rename and move .gff files
        for file in "$folder"/*.gff; do
            if [ -f "$file" ]; then
                mv "$file" "$DEST_DIR/${basename}.gff"
            fi
        done
    fi
done

echo "Files renamed and moved to $DEST_DIR."


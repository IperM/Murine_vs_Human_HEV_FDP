#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_folder output_file.fasta"
    exit 1
fi

input_folder="$1"
output_file="$2"

# Loop through each file in the input folder
for file in "$input_folder"/*; do
    # Check if file exists and is a regular file
    if [ -f "$file" ]; then
        # Search for the desired block
        block=$(sed -n '/>>Overrepresented sequences\twarn/,/^>>/p' "$file")

        # Extract the first and last fields
        first_field=$(echo "$block" | awk 'NR==3 {print $1}')
        last_field=$(echo "$block" | awk '/^>>/ {exit} END {print $NF}')

        # Write to output file in fasta format
        echo ">$first_field" >> "$output_file"
        echo "$last_field" >> "$output_file"
    fi
done

echo "Extraction complete. Output written to $output_file"

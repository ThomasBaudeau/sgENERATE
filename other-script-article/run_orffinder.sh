#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <multi-FASTA file>"
    exit 1
fi

# Check if the file exists
if [ ! -f "$1" ]; then
    echo "Error: File '$1' not found!"
    exit 1
fi

# Loop through the multi-FASTA file
while IFS= read -r line; do
    # Check for the header line
    if [[ $line == ">"* ]]; then
        # Extract the sequence name
        seq_name=$(echo "$line" | cut -d ">" -f 2)
        echo "Extracting sequence: $seq_name"
        # Create a new file for the sequence
        seq_file="${seq_name// /_}.fa"
        echo "$line" > "$seq_file"
    else
        # Append sequence lines to the current sequence file
        echo "$line" >> "$seq_file";
        snakemake -c 1 --use-conda "$seq_name"_COV-ref.exonerate_affine_global.txt;
        snakemake -c 1 --use-conda "$seq_name"___COV-ref.exonerate_est2_genome.txt;
        snakemake -c 1 --use-conda "$seq_name"_COV-ref.exonerate_align_local_S_no.txt;
        rm $seq_file
    fi
done < "$1"

echo "Extraction complete!"


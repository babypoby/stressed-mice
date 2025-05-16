#!/bin/bash
module load stack/2024-06 gcc/12.2.0 bedtools2/2.31.0
# Define the paths
INPUT_DIR="$SCRATCH/data_oxidation/oxidation_per_site"
OUTPUT_DIR="$SCRATCH/data_oxidation/gene_bodies_intersect"
ANNOTATION_FILE="/nfs/nas12.ethz.ch/fs1201/green_groups_let_public/Euler/Vakil/mouse_genome_annotation/Genes_Promoters_CpG_islands_for_Tae/knownGenes_canononical_GRCm39_GENCODE.VM36.bed"

# Check if INPUT_DIR exists
if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist."
    exit 1
fi

# Create the output directory
mkdir -p "$OUTPUT_DIR"

# List the directory contents for debugging
echo "Listing files in $INPUT_DIR:"
ls -la "$INPUT_DIR"

# Process all files that end with .bedgraph
for bedgraph in $(ls "$INPUT_DIR"/*.bedgraph 2>/dev/null); do
    # Extract the base filename without extension
    base_filename=$(basename "$bedgraph" .bedgraph)
    
    # Run bedtools intersect
    echo "Processing $bedgraph..."
    bedtools intersect -a "$bedgraph" \
                     -b "$ANNOTATION_FILE" \
                     -wb -f 1 > "$OUTPUT_DIR/${base_filename}.bed"
    
    # Check if the command was successful
    if [ $? -eq 0 ]; then
        echo "Successfully processed $bedgraph"
    else
        echo "Error processing $bedgraph"
    fi
done

echo "All files have been processed."
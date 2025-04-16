#!/bin/bash

# This script reads the output of sorted files from the previous script and creates an intermediate 
# database of symbolic links to the original files

# Parameters
output_file="sorted_snapshots.txt"
SYM_DIR="./1_SYMBOLIC"

# Read the output file and extract the INIT and END counters
START_COUNTER=$(head -n 1 "$output_file" | awk '{print $1}')
END_COUNTER=$(tail -n 1 "$output_file" | awk '{print $1}')

# Create the directory structure ./1_SYMBOLIC/INIT_END
mkdir -p "$SYM_DIR/${START_COUNTER}_${END_COUNTER}"

# Read the output file line by line and create symbolic links
while read -r line; do
    
    # Extract counter, time, and file path from each line
    counter=$(echo "$line" | awk '{print $1}')
    file_path=$(echo "$line" | awk '{print $3}')

    echo "File $counter"

    # Create the symbolic link in the appropriate folder
    ln -s "../../$file_path" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/la2naca_wing0.f${counter}"

done < "$output_file"

mv "$output_file" "$SYM_DIR/${output_file%.txt}.${START_COUNTER}_${END_COUNTER}.txt"

echo "Symbolic links created in: $SYM_DIR/${START_COUNTER}_${END_COUNTER}/"

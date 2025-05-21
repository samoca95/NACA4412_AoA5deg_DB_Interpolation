#!/bin/bash

# This script reads the output of sorted files from the previous script and creates an intermediate 
# database of symbolic links to the original files

# Parameters
output_file="validate_raw.txt"
SYM_DIR="./2_INTERPOLATED_RAW_SS"

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

    directory_path=$(dirname -- "$file_path")
    file_name=$(basename -- "$file_path")
    file_extension="${file_name##*.}"

    echo "File $counter"

    # Create the symbolic link in the appropriate folder
    ln -s "../../$directory_path/intHistory.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intHistory.f${counter}"
    ln -s "../../$directory_path/intVX.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intVX.f${counter}"
    ln -s "../../$directory_path/intVY.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intVY.f${counter}"
    ln -s "../../$directory_path/intVZ.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intVZ.f${counter}"
    ln -s "../../$directory_path/intOX.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intOX.f${counter}"
    ln -s "../../$directory_path/intOY.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intOY.f${counter}"
    ln -s "../../$directory_path/intOZ.$file_extension" "$SYM_DIR/${START_COUNTER}_${END_COUNTER}/intOZ.f${counter}"

done < "$output_file"

mv "$output_file" "$SYM_DIR/${output_file%.txt}.${START_COUNTER}_${END_COUNTER}.txt"

echo "Symbolic links created in: $SYM_DIR/${START_COUNTER}_${END_COUNTER}/"

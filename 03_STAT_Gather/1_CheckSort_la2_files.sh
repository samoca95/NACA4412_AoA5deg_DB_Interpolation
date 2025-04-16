#!/bin/bash

# This script reads the raw db folders and creates a list of the files ordered by time and 
# where the time separation between them has been validated

# Parameters
INIT_SNAPSHOT_NUM=1  # Starting snapshot number
RAW_DIR="./0_RAW/*"
FILE="stsnaca_wing0.f*"
output_file="sorted_stats.txt"

# Temporary file to store unsorted snd filtered data
tmp_file=$(mktemp)
tmp_filtered=$(mktemp)

# Loop over all matching files in the subdirectories
find -L $RAW_DIR -type f -name $FILE | while read -r filepath; do
    echo "Reading $filepath"    

    # Extract the first 80 chars of the snapshot with the metadata
    header=$(head -c 80 "$filepath")

    # Extract the 8th whitespace-separated field
    time=$(echo "$header" | awk '{print $8}')
    # echo $time

    # Extract res_* number from the path (res_# or res_##)
    if [[ "$filepath" =~ res_([0-9]+)/ ]]; then
        resnum=${BASH_REMATCH[1]}
    else
        echo "WARNING: Could not extract res_* number from $filepath" >&2
        continue
    fi

    # Validate time
    if [[ "$time" =~ ^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
        echo "$time $resnum $filepath" >> "$tmp_file"
    else
        echo "WARNING: Invalid time in $filepath" >&2
    fi
done

# Sort by time then by res number descending (so highest comes first)
sort -k1,1n -k2,2nr "$tmp_file" | \
    awk '!seen[$1]++ {print $1, $3}' > "$tmp_filtered"

# Sort again by time (in case time order was broken by the res selection)
sort -n -k1,1n "$tmp_filtered" > "${tmp_filtered}_sorted"

# # Validate time step uniformity on sorted time list
# awk -v dt="$DT_EXPECTED" '
    # NR == 1 { prev = $1; next }
    # {
        # delta = $1 - prev
        # if ((delta - dt > 1e-7) || (delta - dt < -1e-7)) {
            # printf "ERROR: Non-uniform timestep between %.12f and %.12f (Δt = %.12f, expected %.12f)\n", prev, $1, delta, dt > "/dev/stderr"
            # exit 1
        # }
        # prev = $1
    # }
# ' "${tmp_filtered}_sorted"

# Write result to final output file with snapshot counter
counter=$INIT_SNAPSHOT_NUM
{
    # Add the snapshot counter and file path
    awk -v counter="$counter" '{ printf "%05d %.12f %s\n", counter++, $1, $2 }' "${tmp_filtered}_sorted"
} > "$output_file"

# Clean up
rm "$tmp_file" "$tmp_filtered" "${tmp_filtered}_sorted"

echo "Filtered and validated list saved to: $output_file"

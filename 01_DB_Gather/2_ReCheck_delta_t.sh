#!/bin/bash

filename="sorted_snapshots.txt"
expected_dt=0.006
prev_time=""
line_num=0
errors=0

# Read the file line by line
while IFS= read -r line; do
    ((line_num++))
    
    # Extract the time (second field)
    current_time=$(echo "$line" | awk '{print $2}')
    
    # Extract current snapshot number
    i_field=$(echo "$line" | awk '{print $1}')

    if [ -n "$prev_time" ]; then
        # Calculate the difference using awk for floating point precision
        local_diff=$(awk -v curr="$current_time" -v prev="$prev_time" 'BEGIN {printf "%.9f", curr - prev}')    

        # Compare with expected difference (using awk for floating point comparison)
        if ! awk -v local_diff="$local_diff" -v expected_dt="$expected_dt" 'BEGIN {exit (local_diff != expected_dt)}'; then
            echo "Error at line $line_num: Time difference is $local_diff (expected $expected_dt) between lines $((line_num-1)) and $line_num"
            ((errors++))
        else 
            echo "Field $i_field OK"
        fi
    fi
    
    prev_time=$current_time
done < "$filename"

if [ $errors -eq 0 ]; then
    echo "All time differences are $expected_dt"
else
    echo "Found $errors time difference errors"
fi

exit $errors

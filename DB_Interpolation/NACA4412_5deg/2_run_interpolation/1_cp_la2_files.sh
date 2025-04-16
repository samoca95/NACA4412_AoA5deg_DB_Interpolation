#!/bin/bash

# RUN_FOLDER=$(ls | grep 'Run')


RUN_FOLDER=('/home/oner/oner751450/scratch/NACA4412_5deg_DB/1_SYMBOLIC/04001_05000') #'Run17_11325_12075' 'Run18_12075_12825') 


i_file=1
for i_folder in "${RUN_FOLDER[@]}"; do
    echo "coping files from $i_folder"
    
    for orig_file in "$i_folder"/*; do
        if [[ -f "$orig_file" ]]; then
            new_file=la2naca_wing0.f$(printf "%05d\n" $i_file)
            ln -s $orig_file ./ZSTRUCT/$new_file
            echo "$orig_file -> $new_file"
            i_file=$(($i_file+1))   
        fi
    done
done

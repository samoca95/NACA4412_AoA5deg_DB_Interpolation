#!/bin/bash

# RUN_FOLDER=$(ls | grep 'Run')


RUN_FOLDER=('Run16_10575_11325' 'Run17_11325_12075' 'Run18_12075_12825' 'Run19_12852_13575' 'Run20_13575_14325' 'Run21_14325_15075' 'Run22_15075_15825' 'Run23_15825_16575' 'Run24_16575_17325' 'Run25_17325_18075' 'Run26_18075_18825')


i_sts_file=0

for i_folder in "${RUN_FOLDER[@]}"
do
    echo "cp sts files from $i_folder"

    STATFILE=$i_folder/ptsnaca_wing0.f00001

    if test -f "$STATFILE"; then
      i_sts_file=$(($i_sts_file+1))
      file_name=ptsnaca_wing0.f$(printf "%05d\n" $i_sts_file)
      echo $file_name
	    cp $STATFILE ../../time_series/naca200k_ref/$file_name
    fi

    STATFILE=$i_folder/ptsnaca_wing0.f00002

    if test -f "$STATFILE"; then
      i_sts_file=$(($i_sts_file+1))
      file_name=ptsnaca_wing0.f$(printf "%05d\n" $i_sts_file)
      echo $file_name
      cp $STATFILE ../../time_series/naca200k_ref/$file_name
    fi

done

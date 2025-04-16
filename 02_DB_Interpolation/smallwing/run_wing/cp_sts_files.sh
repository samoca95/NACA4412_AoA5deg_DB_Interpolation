#!/bin/bash

RUN_FOLDER=$(ls | grep 'Run')

i_sts_file=0

for i_folder in $RUN_FOLDER
do
    echo "cp sts files from $i_folder"

    STATFILE=$i_folder/stsnaca_wing0.f00001

    if test -f "$STATFILE"; then
      i_sts_file=$(($i_sts_file+1))
      file_name=stsnaca_wing0.f$(printf "%05d\n" $i_sts_file)
      echo $file_name
	    cp $STATFILE $file_name
    fi

    STATFILE=$i_folder/stsnaca_wing0.f00002

    if test -f "$STATFILE"; then
      i_sts_file=$(($i_sts_file+1))
      file_name=stsnaca_wing0.f$(printf "%05d\n" $i_sts_file)
      echo $file_name
      cp $STATFILE $file_name
    fi

done

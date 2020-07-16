#!/bin/bash

if [ -z "$2" ]
  then
    echo "Please provide sequence and reducing_end_only_flag"
    exit
fi

./submit/submit $1 $2 10980
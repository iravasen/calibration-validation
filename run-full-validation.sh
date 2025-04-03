#!/bin/bash

# Run it like this
# ./run-full-validation.sh <path where the rawTF folders are stored>
BASE_PATH=$1

echo "Validation started: $(date +"%Y-%m-%d %H:%M:%S")"

./run-single-validation.sh 560692 thr_scan ${BASE_PATH}
./run-single-validation.sh 560721 dig_scan ${BASE_PATH}
./run-single-validation.sh 560724 vresetd_1d ${BASE_PATH}
./run-single-validation.sh 560725 vresetd_2d ${BASE_PATH}
./run-single-validation.sh 561718 ps_1d ${BASE_PATH}
./run-single-validation.sh 560727 ps_2d ${BASE_PATH}
./run-single-validation.sh 560728 tot_full ${BASE_PATH}
./run-single-validation.sh 560729 vcasn ${BASE_PATH}
./run-single-validation.sh 560730 ithr ${BASE_PATH}

echo "Validation finished: $(date +"%Y-%m-%d %H:%M:%S")"

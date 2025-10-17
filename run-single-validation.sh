#!/bin/bash 

# Run it like this
# ./run-single-validation.sh <run_number> <scan_name> <base_path_rawTF>
# <scan_name> can be: thr_scan, dig_scan, vresetd_1d, vresetd_2d, ps_1d, ps_2d, tot_full, vcasn, ithr

#####
# Prepare file with list of RAW TF to be processed
#####
BASE_PATH=$3
RUN_NUMBER=$1
SCAN_NAME=$2
find . ${BASE_PATH}/${RUN_NUMBER}_${SCAN_NAME}/ | grep "o2_raw" > ${RUN_NUMBER}.dat

#####
# Run O2 workflows based on the run type and save everything in log file run_number.log
#####

NINJ=50
ADD_OPTIONS=""
if [[ ${SCAN_NAME} == "vresetd_2d" || ${SCAN_NAME} == "ps_2d" ]]; then
	NINJ=10
fi
if [[ ${SCAN_NAME} == "tot_full" ]]; then
	ADD_OPTIONS="--calculate-slope --charge-a 30 --charge-b 60"
	NINJ=10
fi

echo -e "Running workflow on run ${RUN_NUMBER} - ${SCAN_NAME}"

o2-raw-tf-reader-workflow -b --input-data "${RUN_NUMBER}.dat" --onlyDet ITS --shm-segment-size 30000000000 | \
o2-itsmft-stf-decoder-workflow -b --rof-lenght-error-freq -1 --always-parse-trigger --allow-empty-rofs --decoder-verbosity 0 --ignore-noise-map --nthreads 1 --pipeline its-stf-decoder:2 --enable-calib-data --digits --no-clusters --no-cluster-patterns | \
o2-its-threshold-calib-workflow -b --s-curve-col-step 256 --ninj ${NINJ} ${ADD_OPTIONS} --ccdb-mgr-url "http://alice-ccdb.cern.ch" --enable-single-pix-tag --chip-mod-selector 0 --chip-mod-base 3 --meta-output-dir "./" --nthreads 1 --fittype derivative | \
o2-its-threshold-calib-workflow -b --s-curve-col-step 256 --ninj ${NINJ} ${ADD_OPTIONS} --ccdb-mgr-url "http://alice-ccdb.cern.ch" --enable-single-pix-tag --chip-mod-selector 1 --chip-mod-base 3 --meta-output-dir "./" --nthreads 1 --fittype derivative | \
o2-its-threshold-calib-workflow -b --s-curve-col-step 256 --ninj ${NINJ} ${ADD_OPTIONS} --ccdb-mgr-url "http://alice-ccdb.cern.ch" --enable-single-pix-tag --chip-mod-selector 2 --chip-mod-base 3 --meta-output-dir "./" --nthreads 1 --fittype derivative | \
o2-its-threshold-aggregator-workflow -b > ${RUN_NUMBER}.log

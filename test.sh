#!/bin/sh
./hdr.sh configs/configAD.txt
./hdr.sh configs/configAD_with_region_file.txt
./hdr.sh configs/configAD_with_include_vcftable.txt
./hdr.sh configs/configAD_with_exclude_vcftable.txt
./hdr.sh configs/configAD_with_include_exclude_vcftable.txt
./hdr.sh configs/configAD_run_with_vcftable_from_configAD.txt
./hdr.sh configs/configAR.txt
./hdr.sh configs/configAR_with_big_window.txt
./hdr.sh configs/configDD.txt
./hdr.sh configs/configDD_run_genomecov.txt
./hdr.sh configs/configDD_run_select_minmaxregion_vcftable.txt
./hdr.sh configs/configDD_run_both_findrun_and_genomecov.txt

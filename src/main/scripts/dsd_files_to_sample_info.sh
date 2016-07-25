#!/bin/bash
#
# Wrapper script to run DSD sample meta data and pedigree script
#

function load_config() {
    BASE=../..
    CONFIG=`sed 's/\/\/.*$//' $BASE/pipeline/config.groovy` 
    eval "$CONFIG"
}

load_config

if [ ! -e data ] || [ ! -e ../../batches ];
then
	echo
	echo "Please run this script from within the batch directory for a DSD Haloplex run, after copying FASTQ files to the data directory"
	echo
	exit 1
fi

$GROOVY -cp ../../tools/groovy-ngs-utils/1.0.2/groovy-ngs-utils.jar  ../../tools/haloplex_dsd/src/main/scripts/dsd_files_to_sample_info.groovy $*

#!/bin/bash

# Note: this script is meant to serve as a repository of the commands
# to use for processing barcode sequencing data. It is not meant to
# be run all at once. Please run the steps individually and make
# sure each completes successfully before the next step is run.

# Set the location of this barseq_counter folder
export BARSEQ_DIR=/project/csbio/chemical_genomics_pipeline/barseq_counter_minimal

# Set path to compiled agrep executable
# We have included a linux-compiled version in the "agrep-3.41" folder
export AGREP=$BARSEQ_DIR/agrep/agrep

# Prepend the barseq_counter scripts directory to your path, so you can
# call the commands by just their names.
# NOTE: this will not work if the $BARSEQ_DIR path does not exist (i.e.
# you typed it in wrong).
pathPrepend() {
    if [ -d "$1" ] && [[ ":$PATH:" != *":$1:"* ]]; then
        PATH="$1:${PATH:+"$PATH"}"
    fi
}

pathPrepend $BARSEQ_DIR/scripts
export PATH

# Run the scripts
# CD into your working directory
# In this case, we'll just use the example folder
cd $BARSEQ_DIR/example

# Download the example data and deposit in the "data" folder
# Please allow some time for this to download
mkdir -p data
wget http://lovelace-umh.cs.umn.edu/chemical_genomics_tools/barseq_counter/sample_dataset/L1.fastq ./data
wget http://lovelace-umh.cs.umn.edu/chemical_genomics_tools/barseq_counter/sample_dataset/L2.fastq ./data

# This step grabs the index tags and barcodes from sequences with the common priming site
preprocess_MIseq_10bp.py data barseq.txt

# This step maps barcodes and index tags to a strain by condition count matrix 
processBarSeq_rd.py barseq.txt decode.txt $BARSEQ_DIR/barcodes/allupbarcodes.txt barseq.processed.cerevisiae.txt

# Create a count matrix with yeast common names
convertORF2common.py barseq.processed.cerevisiae.txt barseq.processed.cerevisiae.common.txt

# This step generates reports
mkdir -p example_reports
generateReport_MIseq.py data barseq.txt barseq.processed.cerevisiae.txt example_reports/cerevisiae

# Then, perform statistical analysis using EdgeR
# Note: if the user's Rscript is not located at "/usr/bin/Rscript",
# then they will need to change the location in the runEdgeR.R script.
runEdgeR.R --threshold 10 barseq.processed.cerevisiae.txt controls.txt


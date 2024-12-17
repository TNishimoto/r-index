#!/bin/bash
cd `dirname $0`


SCRIPT_DIR=$(cd $(dirname $0); pwd)
echo "Script directory: $SCRIPT_DIR"

#INPUTDIR=../data/external/repcorpus
#OUTPUTDIR=../data/bwt

INPUTDIR=/Users/nishimoto/Documents/test_data/bwt
OUTPUTDIR=$SCRIPT_DIR/../build/indexes


for file in $INPUTDIR/*; do
    if [ -f "$file" ]; then
        echo "Full path: $(realpath "$file")"
        echo "File name: $(basename "$file")"
      $SCRIPT_DIR/../build/mod_ri-build -o $OUTPUTDIR/$(basename "$file") $INPUTDIR/$(basename "$file")
    fi
done

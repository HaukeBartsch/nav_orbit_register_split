#!/usr/bin/env bash
# 
# Call 3DSlicer with a volume and a label file
# (based on training in dataset.json).
#
# use 0 as the default
pair="0"

if [ $# -eq 1 ]; then
    pair="$1"
fi

# open slicer and load one example
image=$(jq -r ".training[$pair].image" metadata.json)
mask=$(jq -r ".training[$pair].label" metadata.json)
image=$(readlink -f ${image})
mask=$(readlink -f ${mask})
image_name=$(basename "$image")
mask_name=$(basename "$mask")

/Applications/Slicer.app/Contents/MacOS/Slicer --python-code "slicer.util.loadVolume('${image}',{'name': '${image_name}'}); slicer.util.loadSegmentation('${mask}',{'name': '${mask_name}'})"
#
# To safe the output back convert to a binary label first.
#

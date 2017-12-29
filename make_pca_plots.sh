#!/bin/sh

cd preprocess
python extract_diffs.py
python gen_oneline_format.py
cd ..

cd count
python count_derived_each_lineage_nosingle.py
cd ..

cd plot
python make_finescale_PCA_nosingle.py
cd ..

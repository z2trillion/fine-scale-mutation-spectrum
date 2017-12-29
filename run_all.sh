#!/bin/sh

cd preprocess
python extract_diffs.py
python gen_oneline_format.py
cd ..

cd count
python process_chromosomes.py;
cd ..

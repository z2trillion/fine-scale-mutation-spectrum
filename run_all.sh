#!/bin/sh

cd preprocess
python extract_diffs.py
python gen_oneline_format.py
cd ..

cd count
python get_finescale_mut_spectra_nosingle.py 22;
python get_finescale_phyloP_mut_type_v_freq_nosingle.py 22;
python get_finescale_repeats_mut_type_v_freq_nosingle.py 22;
cd ..

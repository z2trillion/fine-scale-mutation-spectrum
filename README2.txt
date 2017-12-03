To generated mutation spectrum summary files from the 1000 Genomes VCFs, run:
   python get_finescale_mut_spectra_nosingle.py chrom
   python get_finescale_repeats_mut_type_v_freq_nosingle.py chrom
   python get_finescale_phyloP_mut_type_v_freq_nosingle.py chrom
for every autosome number "chrom" ranging from 1 to 22

This will populate an empty directory "finescale_mut_spectra/" with all of the summary files needed to generate a heat map comparison between East Asian populations. For reference, these summary files are already included. You can produce a heatmap visualization using:
     python make_heatmap_ASN_filtered.py 

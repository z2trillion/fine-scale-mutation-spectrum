Before you generate mutation spectrum summary files from the 1000 Genomes VCFs, you will need to put the human reference sequence and the human/chimp alignment into a certain format.

To reformat the human reference sequence, copy all fasta files "chrN.fa", N ranging from 1 to 22, into the folder "hg19_reference". For each N, run:
   python gen_oneline_format.py N 
to generate the file "chrN_oneline.txt". For reference, the file "chr22_oneline.txt" is already included in the hg19_reference folder.

Similarly, you will need to extract a list of human/chimp differences from each file "chrN.hg19.panTro.net.axt.gz". To do so, copy the chain files into the folder "hg19_chimp_align" and run
	   python extract_diffs.py N
for each chromsome number N to generate human_chimp_diffs_chrN.txt. The output file "human_chimp_diffs_chr22.txt" is included for reference.

Once the reference and ancestral state files are made, to generated mutation spectrum summary files from the 1000 Genomes VCFs, run:
   python get_finescale_mut_spectra_nosingle.py N
   python get_finescale_repeats_mut_type_v_freq_nosingle.py N
   python get_finescale_phyloP_mut_type_v_freq_nosingle.py N
for every autosome number "N" ranging from 1 to 22.

I would love the package that replaces this pipeline to be able to input any number of BED files and exclude the regions masked by those BED files from the mutation spectrum analysis. At the moment, my pipeline excludes 2 files: conserved regions and annotated repeats.

This will populate an empty directory "finescale_mut_spectra/" with all of the summary files needed to generate a heat map comparison between East Asian populations. For reference, the output you should generate is already included in the folder "finescale_mut_spectra_KH." You can produce a heatmap visualization of the mutation spectrum differences between populations using:
     python make_heatmap_ASN_filtered.py 

(If you want to make the heat map from the data I preprocessed instead of new files that you generate yourself, replace the folder name "finescale_mut_spectra" with the folder name "finescale_mut_spectra_KH”)

To make the PCA plot in Figure 4, you first run 

python count_derived_each_lineage_nosingle.py N

for each chromosome number "N" to produce the files "derived_each_lineage_chrN_nosingle.txt." This depends on the hg19 and chimp alignment files you already generated to make the heatmaps, as well as the 1000 Genomes VCFs. Once you have these files, you can run

python make_finescale_PCA_nosingle.py

which will generate a PCA for each 1000 Genomes population. 

To make the other PCA plot in Figure 1, which combines all the data, run:

python make_PCA_diploid_nosingle.py

As you may have noticed, the human data contain continent-scale labels and more specific population-scale labels. For the general package, I’d like the user to be able to input any population labeling of their samples they want to look at. 
# Visualizing the footprints of mutation spectrum evolution
These are the lecture notes from Kelley Harris' talk at the 2018 Workshop on
Population and Speciation Genomics in Český Krumlov.

## Background
In today’s lecture, you learned that the mutational process is not completely
clock-like, but is better regarded as a mixture of many different processes that
damage DNA and cause it to be replicated unfaithfully. These different processes
can be disentangled by looking at the relative rates of different types of
mutations, and comparing these relative rates across populations or species. For
more background, see
[Harris & Pritchard eLife 2017](https://elifesciences.org/articles/24284).  We
will be using Python scripts to extract mutation spectrum information from a VCF
of human variation data, the final 1000 Genomes Phase 3 panel, using the human
hg19 reference genome sequence to identify the sequence context in which each
variant occurs and using a human-chimpanzee genome alignment for ancestral
allele polarization. We will then visualize these differences in a two different
ways: using a heat map and a principal components analysis plot.  If you have
time at the end of the exercise or after the workshop, consider modifying the
pipeline to have a look at mutation spectrum variation in a whole genome dataset
of your own. Very little is known about the ancestral dynamics of the mutation
process in species outside of humans and great apes, and it would be exciting to
see how mutagenesis has been varying in species with different demographic
histories and reproductive strategies.

## Getting started
The python scripts that you will need to do the exercise are provided in a
folder called “fine-scale-mutation-spectrum-master.” You will notice that the
folder contains a directory called “data,” which contains processed data for
most of the human genome but is missing chromosome 22. Your first step will be
to download the missing chromosome 22 data from the UCSC genome browser and the
1000 Genomes FTP server, so that you get to practice looking through the rich
treasure troves of public data available on the internet.

VCF files list the SNPs that occur in a set of whole genome sequences. For
each SNP, they list the position in the reference sequence where the occurs,
the reference base at that position, and the alternate base that occurs in
many of the reads that map to that position. However, the VCF file is missing
two pieces of information that are crucial for mutation spectrum analysis. One
piece is ancestral state information: is the reference allele or the alternate
allele the one that more recently arose by mutation? The other missing piece
is the triplet context: which base pairs flank the mutation in the reference
sequence? To add these pieces of information to the VCF, we must download but
the reference genome sequence and an alignment between the reference genome
and an outgroup genome, then reformat these files.

## Downloading and preprocessing your reference files
Your first task is to use the UCSC Genome Browser to find the hg19 human
reference sequence and an alignment between hg19 and the most up-to-date
chimpanzee reference sequence. Download the reference sequence and the
hg19-chimp reference sequence for chromosome 22. Put these sequences into the
directories “data/hg19_reference” and data/hg19_chimp_align”, respectively.

We will also use the UCSC Genome Browser to obtain bed file annotations of
regions of the genome that we may want to filter out of our analysis. We want to
filter out the regions annotated as “conserved” by the program PhastCons, which
appear to be evolving under significant amounts of natural selection, and we
also want to filter out regions annotated as repeats by RepeatMasker. First,
discuss with your neighbor why it might be a good idea to filter out these
regions. Next, use the Genome Browser website to look for the files
“phastConsElements100way.txt” and “nestedRepeats.txt”, making sure to get the
files that are indexed with respect to the hg19 reference. Put these files into
the empty folder “data/bed”.

The directory “preprocess” contains a script to reformat the alignment to the
outgroup. You can run it from that directory with this command:
```
python extract_diffs.py [chromsome_number]
```

## Extracting the mutation spectrum from the VCF
Armed with your preprocessed reference files, you can now extract the mutation
spectrum of chromosome 22. I pre-extracted the mutation spectra of the other
chromosomes; you can see that these processed spectra are located in the folder
“finescale_mut_spectra”, and that the chromosome 22 mutation spectrum is
missing.

To generate the chromosome 22 mutation spectrum, you will need to google the
1000 Genomes Project Phase 3 FTP server, download the VCF file for chromosome
22, and add it to the empty folder “data/vcfs”.  (Hint: look for the 1000
Genomes FTP, look through the directory structure for a folder called “release,”
then look inside the release folder for the subdirectory that has been modified
most recently. This will contain the VCF you want). Then, use the following
command from the “count” directory to count the number of SNPs occurring in each
triplet sequence context at each allele frequency:
```
python process_chromosomes.py [chromosome_number]
```
Note that singletons are excluded from the counts due to concerns that they
could be enriched for sequencing errors.

You should find that the directory “finescale_mut_spectra” now contains new with
data from chromosome 22. Some tabulate all mutations that occur in chromosome
22, for each population separately, while the others tabulate mutations that
occur only in conserved regions or only within repeat elements. When you start
plotting mutation spectrum data, the plotting scripts will subtract away the
mutations that occur within the regions you want to filter away.

It is easy to modify the process_chromosomes.py script to output mutation
spectra contained within any genomic annotation you can specify with a bed file,
say, with high recombination rate regions or within open chromatin. If you have
time at the end of the exercise, you can experiment with doing this.  

## Creating a mutation spectrum heat map
Now that you’ve tabulated mutation spectrum data for each population, it’s time
to visualize the data so we can pick out differences between populations. The
script you will use is “make_heatmap.py”, within the “plot” directory. There are
a lot of different options built in, which I hope you’ll play around with.

First, make a figure with three side by side mutation spectrum comparisons:
Great Britain (GBR) vs Chinese from Beijing (CHB); Finnish (FIN) vs Japanese
from Tokyo (JPT); GBR vs FIN; JPT vs CHB. The command for doing this is:
```
python make_heatmap.py -p GBR CHB FIN JPT GBR FIN
```
Each grid square represents a mutation type (a mutation from some triplet ABC to
some other triplet ADC). A red square means that the population listed first has
a higher proportion of ABC->ADC mutations, relative to other mutation types,
than the population listed second. A blue square means that the population
listed first is enriched for ABC->ADC mutations relative to the population
listed second. In both cases, a white dot on the square means that the
difference between populations is significant (chi-square p-value < 1e-5).

Take a moment to compare your heat map with your neighbor’s and discuss your
observations about it.

Once you’ve successfully your first heat map and discussed it, it’s time to
start playing with the following optional command line arguments:
```
-c chrom1 chrom2 ... : only use mutation counts from the chromosomes listed.
                       Default is range(1,23).
-f min_frequency max_frequency: Only count mutations whose frequencies range
                                between the listed min and max, inclusive.
                                Default is min_frequency=0, max_frequency=1.
-e : If this flag is used, exclude the SNPs occurring in PhyloP conserved regions
     and repetitive regions.  
-i : Analyze chromosomes individually: If this flag is used, output a separate
     heat map for each specified chromosome individually.
-p pval : Print white dots on all squares representing mutation types that are
          significantly enriched in one population at a chi-square significance
          level of pval. Default = 1e-5.
```
Use these flags to generate the following heat maps:
- Generate a heat map using only chromosome 22 data with the following 4
comparisons side by side: Great Britain vs Han Chinese from Beijing; Finland vs
Japanese; Great Britain vs Finland; Japanese vs Han Chinese from Beijing.
Discuss.
- Try generating the same heat map without excluding repeats or conserved regions.
Do you see much of a difference?
- Repeat using precomputed summary data for all 22 autosomes.
- Generate another heat map with the same population comparisons (all 22
autosomes), but restricting to alleles with frequency less than 0.98. Discuss
the differences.
- Generate two more heat maps with the same population comparisons (all 22
autosomes), one restricting to alleles with frequency 0.5-0.9; one with
frequencies 0.1-0.5; one with frequency <0.1. Discuss the differences.
- Generate a new heat map plot with 4 panels that each compare Great Britain vs
Han Chinese from Beijing, each panel using data from a different chromosome. The
panels should correspond to chromosomes 1, 2, 22, and X. What do you observe?

## Principal component analysis
Although the mutation spectrum heat map is good for pinpointing differences
between populations in the types of mutations that have been accumulating, it is
less suitable for visualizing how the differentiation between populations
compares to the heterogeneity within populations. To explore how much the
mutation spectrum varies within populations, we will use a type of principal
component analysis (PCA).

You have probably seen population-differentiation PCAs before, as pioneered in
[Novembre, et al. Nature 2008](https://www.nature.com/articles/nature07331).
These PCAs are generated by summarizing each genome as a 10,000(+)-dimensional
genotype vector, then projecting these vectors onto the optimal 2-dimensional
plane for viewing. Here, we will summarize each genome as a 96-dimensional
vector that counts the proportion of derived alleles from that genome falling
into the 96 possible frequency classes (You might think that there are 192
mutation categories of the form ABC->ADC, but recall that each mutation type is
the same as its reverse complement, e.g. CCA->CTA is the same as TGG->TAG.)

To make mutation spectrum PCA plots, we need to count mutations in a different
way than we counted them for the heat maps. Instead of counting the number of
mutations that occur in each frequency class, we need to count the number that
occur within each genome. As before, we exclude singletons due to data quality
concerns. To count the mutation types in each genome for chromosome 22, run the
following command in the “count” directory:
```
python count_derived_each_lineage_nosingle.py -c 22
```
As before, counts for the remaining chromosomes are waiting precomputed for you
in the “finescale_mut_spectra” folder.

After generating mutation counts for each genome, we need to switch over to the
“plot” directory and start plotting PCAs. The command
```
python make_finescale_PCA_nosingle.py YRI GBR
```
will make a PCA showing how distinct the Yoruban from Ibadan (YRI) genomes are
from the Great Britain (GBR) genomes in terms of mutation spectrum, compared to
the variation within the two groups. As with heat maps, you can restrict to a
subset of chromosomes using the -c flag:
```
-c chrom1 chrom2 ... : only use mutation counts from the chromosomes listed.
                       Default is range(1,23).
```
Here are some exercises to get you acquainted with mutation spectrum PCAs:
- Using the chromosome 22 data you just generated, make a PCA plot comparing the
following groups: Yoruba (YRI), African Americans from the Southwest US (ASW),
Mexicans from LA (MXL), Han Chinese from Shanghai (CHS), Japanese (JPT), Great
Britain (GBR), Finnish (FIN).
- Make another PCA for the same groups using all 22 autosomes. Why do these look
different?
- Generate continent-specific PCAs: one for the set of all East Asian groups (JPT,
CHS, CHB, CDX, KHV) and a second PCA for the set of all Europeans (CEU, GBR,
FIN, TSI, IBS). Discuss the similarities and differences between these plots and
the intercontinental PCAs you generated.
- For your intercontinental and intracontinental PCAs, experiment with including
different numbers of chromosomes and see how much data you need for structure to
emerge.

## Bonus Problem!
If you get tired of exploring patterns of human mutation spectrum variation, you
can choose to accept the challenge of mapping mutation spectrum variation in
Arabidopsis thaliana! To do so, you should first make a copy of the project
folder fine-scale-mutation-spectrum-master. You’ll then need to replace the
human data files in this folder with Arabidopsis data files, as follows:

The Arabidopsis 1001 Genomes Project is a lot like the human 1000 genomes
project in that it sampled individuals from all over the world. In this project,
individuals are called “accessions.” In your copied folder, replace the human
1000 genomes vcf with the vcf of chromosome 3  from Angela Hancock’s exercise

To process the SNPs in this file, you will need to download the A. thaliana
reference chromsome fasta for chromosome 3 from the UCSC browser. Instead of the
alignment between human and chimp, you will need the alignment of A. thaliana
chromosome 3 to its relative A. lyrata, available here:
http://pipeline.lbl.gov/data/araTha04_Araly1

Finally, you will need to replace the file 1000genomes_phase3_sample_IDs.txt
with the accession label file from Angela Hancock’s exercise. Use the “grep”
command to find all instances of 1000genomes_phase3_sample_IDs.txt in the folder
counting scripts and replace it with the name of your new population file.

While you have the labels.py file open, edit the header section that currently
reads:
```
group_to_populations = {
  'EAS': ['CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'CHD'],
  'EUR': ['CEU', 'TSI', 'GBR', 'FIN', 'IBS'],
  'AFR': ['YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ACB', 'ASW'],
  'SAS': ['GIH', 'PJL', 'BEB', 'STU', 'ITU'],
  'AMR': ['CLM', 'MXL', 'PUR', 'PEL'],
}
```
This current header assigns the 3-letter codes that stand for human populations
to the 3-letter codes that stand for larger groups (for example, the first line
says that CHB (Han Chinese from Beijing) and JPT (Japanese from Tokyo) are both
EAS (East Asian) populations. Delete this header and replace it with the
following one, which classifies Arabidopsis populations into the appropriate
continents:
```
group_to_populations = {
  'east_eurasia': ['central_europe', 'central_asia', 'italy_balkan_caucasus'],
  'sweden': ['northern_sweden', 'southern_sweden'],
  'west_eurasia': ['western_europe', 'germany', 'iberian_peninsula'],
  'relict': ['relict'],
  'admixed': ['admixed'],
}
```
Use what you learned from the human mutation spectrum exercise to preprocess the
Arabidopsis data, count mutations, and investigate mutation spectrum variation
across worldwide accessions of Arabidopsis.

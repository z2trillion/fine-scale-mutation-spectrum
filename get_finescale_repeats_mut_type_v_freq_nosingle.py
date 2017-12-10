import sys
from labels import populations
from get_finescale_phyloP_mut_type_v_freq_nosingle import get_finescale
from common import get_conserved

if __name__ == '__main__':
    chrom=sys.argv[1]
    output = {population: 'Ref Alt \n' for population in populations}
    outfile_path = 'finescale_mut_spectra/inrepeats_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
    conserved = get_conserved('data/nestedRepeats.txt', chrom)
    get_finescale(outfile_path, chrom, output, conserved)

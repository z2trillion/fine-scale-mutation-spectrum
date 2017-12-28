import sys
from labels import populations
from get_finescale_phyloP_mut_type_v_freq_nosingle import get_finescale

if __name__ == '__main__':
    chrom=sys.argv[1]
    output={population: 'Mut\n' for population in populations}
    outfile_path = '../finescale_mut_spectra/mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
    conserved = [[0, 1e12]]
    get_finescale(outfile_path, chrom, output, conserved)

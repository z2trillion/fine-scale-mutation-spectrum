import sys
from get_finescale_phyloP_mut_type_v_freq_nosingle import get_finescale

if __name__ == '__main__':
    chrom=sys.argv[1]
    get_finescale('nestedrepeats', chrom)

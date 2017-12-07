import sys
from labels import populations
from get_finescale_phyloP_mut_type_v_freq_nosingle import get_finescale

if __name__ == '__main__':
    chrom=sys.argv[1]
    output = {population: 'Ref Alt \n' for population in populations}
    get_finescale('nestedrepeats', chrom, output)

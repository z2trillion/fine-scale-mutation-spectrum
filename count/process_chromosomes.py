from labels import populations
from common import open_infile, get_conserved, get_chromosomes_from_args
from mutation_counter import MutationCounter, IncludedRegion


def get_finescale(mutation_counter):
    infile, line = open_infile(mutation_counter.chrom)

    mutation_counter.configure(line)

    for line_counter, line in enumerate(infile):
        mutation_counter.process_line(line)

    mutation_counter.write_output()


if __name__ == '__main__':
    chromosomes = get_chromosomes_from_args()
    for chrom in chromosomes:
        chrom = str(chrom)
        included_regions = []

        output={population: 'Mut\n' for population in populations}
        outfile_path = '../finescale_mut_spectra/mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = [[0, 1e12]]
        included_regions.append(IncludedRegion(chrom, output, outfile_path, conserved))

        output = {population: 'Ref Alt \n' for population in populations}
        outfile_path = '../finescale_mut_spectra/inrepeats_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = get_conserved('../data/bed_files/nestedRepeats.txt', chrom)
        included_regions.append(IncludedRegion(chrom, output, outfile_path, conserved))

        output = {population: 'Ref Alt \n' for population in populations}
        outfile_path = '../finescale_mut_spectra/phyloP_conserved_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = get_conserved('../data/bed_files/phastConsElements100way.txt', chrom)
        included_regions.append(IncludedRegion(chrom, output, outfile_path, conserved))

        mutation_counter = MutationCounter(chrom, included_regions)
        get_finescale(mutation_counter)

        print 'finished chrom', chrom

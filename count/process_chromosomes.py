from labels import populations
from common import open_infile, get_conserved, get_chromosomes_from_args
from mutation_counter import MutationCounter


def get_finescale(mutation_counters):
    infile, line = open_infile(mutation_counters[0].chrom)
    for mutation_counter in mutation_counters:
        mutation_counter.configure(line)

    for line_counter, line in enumerate(infile):
        for mutation_counter in mutation_counters:
            mutation_counter.process_line(line)

    for mutation_counter in mutation_counters:
        mutation_counter.write_output()


if __name__ == '__main__':
    chromosomes = get_chromosomes_from_args()
    for chrom in chromosomes:
        chrom = str(chrom)

        output={population: 'Mut\n' for population in populations}
        outfile_path = '../finescale_mut_spectra/mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = [[0, 1e12]]
        everything_mutation_counter = MutationCounter(chrom, output, outfile_path, conserved)

        output = {population: 'Ref Alt \n' for population in populations}
        outfile_path = '../finescale_mut_spectra/inrepeats_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = get_conserved('../data/bed_files/nestedRepeats.txt', chrom)
        nested_repeats_mutation_counter = MutationCounter(chrom, output, outfile_path, conserved)

        output = {population: 'Ref Alt \n' for population in populations}
        outfile_path = '../finescale_mut_spectra/phyloP_conserved_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
        conserved = get_conserved('../data/bed_files/phastConsElements100way.txt', chrom)
        phylop_mutation_counter = MutationCounter(chrom, output, outfile_path, conserved)

        get_finescale([
            everything_mutation_counter,
            nested_repeats_mutation_counter,
            phylop_mutation_counter,
        ])

        print 'finished chrom', chrom

import sys
import argparse
import gzip

from mutations import mutations
from labels import sample_id_to_population, populations


def get_chromosomes_from_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=range(1, 23))

    chromosomes = parser.parse_args(sys.argv[1:]).chromosomes

    valid_chromosomes = ['X'] + [str(i) for i in range(1, 23)]
    for chrom in chromosomes:
        assert chrom in valid_chromosomes
    return chromosomes


def reference_sequence(chromosome_number):
    file_path = '../data/hg19_reference/chr'+chromosome_number+'_oneline.txt'
    with open(file_path) as infile:
        return infile.read()


def get_human_chimp_differences(chromosome_number):
    human_chimp_differences = {}

    file_path = (
        '../data/hg19_chimp_align/human_chimp_diffs_chr' +
        chromosome_number +
        '.txt'
    )
    with open(file_path) as infile:
        infile.next()

        for line in infile:
            if 'SNP' not in line:
                continue
            position, _, _, chimp_allele = line.split()
            human_chimp_differences[int(position)] = chimp_allele

    return human_chimp_differences


def get_column_indices(column_labels):
    population_to_column_indices = {
        population: [] for population in populations
    }

    sample_ids = column_labels[9:]
    for i, sample_id in enumerate(sample_ids):
        population_to_column_indices[
            sample_id_to_population[sample_id]
        ].append(i + 9)

    return population_to_column_indices


def get_column_index_to_population(column_labels):
    sample_ids = column_labels[9:]
    return {
        i + 9: sample_id_to_population[sample_id]
        for i, sample_id in enumerate(sample_ids)
    }


def initialize_mut_count(population_to_column_indices):
    mut_count = {}
    for pop in populations:
        for mut in mutations:
            for i in range(1, 2*len(population_to_column_indices[pop])+1):
                mut_count[(mut, pop, i)] = 0
    return mut_count


def write_output(output, outfile_path, indices, mut_count):
    for pop in populations:
        for mut in mutations:
            output[pop] += mut[0] + '_' + mut[1]
            for i in range(1, 2 * len(indices[pop]) + 1):
                output[pop] += ' ' + str(mut_count[(mut, pop, i)])
            output[pop] += '\n'
        outfile = open(outfile_path % pop, 'w')
        outfile.write(output[pop])
        outfile.close()


def open_infile(chrom):
    file_path = (
        '../data/vcfs/ALL.chr' +
        chrom +
        '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    )
    infile = gzip.open(file_path)

    line = infile.readline()
    while not line.startswith('#CHROM'):
        line = infile.readline()

    return infile, line


def get_conserved(infile_path, chrom):
    infile = open(infile_path)
    lines = infile.readlines()
    infile.close()

    ind = 0
    s = lines[ind].split('\t')
    while not s[1] == 'chr' + chrom:
        ind += 1
        s = lines[ind].split('\t')

    conserved = [(int(s[2]), int(s[3]))]
    ind += 1
    s = lines[ind].split('\t')

    while ind < len(lines) - 1 and s[1] == 'chr' + chrom:
        if int(s[2]) == conserved[-1][-1] + 1:
            new_tup = (conserved[-1][0], int(s[3]))
            conserved.pop()
            conserved.append(new_tup)
        else:
            conserved.append((int(s[2]), int(s[3])))
        ind += 1
        line = lines[ind]
        s = line.strip('\n').split('\t')

    return conserved

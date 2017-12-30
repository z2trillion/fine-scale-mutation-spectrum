import gzip

from itertools import product

from labels import (
    groups,
    population_to_group,
    sample_id_to_population,
)
from mutations import bases
from common import get_chromosomes_from_args


class BadDataQualityError(Exception):
    pass


def parse_line(line):
    (
        chromosome_number,
        position,
        _,  # SNP id
        reference_allele,
        alternate_allele,
        _,  # quality score,
        filter_result,
        info,
        _,  # format (what ever this is???)
        haplotypes
    ) = line.split(None, 9)

    if (
        reference_allele not in list(bases) or
        alternate_allele not in list(bases) or
        filter_result != 'PASS'
    ):
        raise BadDataQualityError

    position = int(position)
    alleles = haplotypes[::2]  # remove '\t' and '|' separators
    derived_count = int(info.split(';')[0].split('=')[-1])

    return reference_allele, alternate_allele, position, derived_count, alleles


def get_mutation(position, reference_allele, alternate_allele, refseq, chimp_alleles):
    context = refseq[position - 2 : position + 1]
    if 'N' in context:
        raise BadDataQualityError

    if chimp_alleles.get(position) == alternate_allele:
        derived_allele='0'
        this_mut=(context[0] + alternate_allele + context[2], reference_allele)
    else:
        derived_allele='1'
        this_mut=(context, alternate_allele)

    return derived_allele, this_mut


def process_line(line, refseq, chimp_alleles):
    (
        reference_allele,
        alternate_allele,
        position,
        derived_count,
        alleles
    ) = parse_line(line)
    derived_allele, this_mut = get_mutation(position, reference_allele,
        alternate_allele, refseq, chimp_alleles)
    return derived_allele, derived_count, this_mut, alleles


def update_counts(alleles, mutation_counts, derived_count, derived_allele, n_lineages, this_mut):
    if derived_count>1 and n_lineages-derived_count>1:
        if derived_allele == '0':
            derived_count=n_lineages-derived_count
        der_observed=0

        for i, allele in enumerate(alleles):
            if allele == derived_allele:
                mutation_counts[(this_mut, i)] += 1
                der_observed += 1
        assert der_observed == derived_count


def write_output(mutation_counts, sample_ids, mutations, n_lineages, chrom):
    output='Mut_type '
    output += ' '.join(
        [sample_id_to_population[sample_id] for sample_id in sample_ids]
    )
    output+='\n'

    for mut in mutations:
        output+=mut[0]+'_'+mut[1]
        for i in range(n_lineages):
            output+=' '+str(mutation_counts[(mut,i)])
        output+='\n'

    with open('../finescale_mut_spectra/derived_each_lineage_chr'+chrom+'_nosingle.txt','w') as outfile:
        outfile.write(output)


def process_chromosome(chrom):
    with open('../data/hg19_reference/chr'+chrom+'_oneline.txt') as infile:
        refseq=infile.read()

    with open('../data/hg19_chimp_align/human_chimp_diffs_chr'+chrom+'.txt') as infile:
        infile.next()
        anc_lines = [line.split() for line in infile if 'SNP' in line]

    chimp_alleles = {}
    for position, _, _, chimp_allele in anc_lines:
        chimp_alleles[int(position)] = chimp_allele

    mutations=[]
    for trimer in product(bases, bases, bases):
        for base in bases:
            if trimer[1] != base:
                mutations.append((''.join(trimer), base))

    gzip_path = '../data/vcfs/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
    with gzip.open(gzip_path) as infile:

        for line in infile:
            if line.startswith('#CHROM'):
                break

        sample_ids = line.split()[9:]
        n_lineages = 2 * len(sample_ids)

        mutation_counts = {
            (mutation, haplotype_index): 0
            for haplotype_index in range(n_lineages)
            for mutation in mutations
        }

        for counter, line in enumerate(infile):
            try:
                derived_allele, derived_count, this_mut, alleles = process_line(line, refseq, chimp_alleles)
            except BadDataQualityError:
                continue

            update_counts(alleles, mutation_counts, derived_count, derived_allele, n_lineages, this_mut)

    write_output(mutation_counts, sample_ids, mutations, n_lineages, chrom)


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=int, nargs='+',
                        default=range(1, 23))

    chromosomes = parser.parse_args(sys.argv[1:]).chromosomes
    for chrom in chromosomes:
        assert 1 <= chrom and chrom <= 22, ('Chromosome %i is unlikely to exist'
                                            % chrom)

    for chrom in chromosomes:
        chromosome = str(chromosome)
        process_chromosome(chromosome)

        print 'finished chrom', chromosome

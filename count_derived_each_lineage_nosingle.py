from copy import deepcopy
import sys
import gzip
import json
from itertools import product
from labels import (
    groups,
    population_to_group,
    sample_id_to_population,
)

chrom = sys.argv[1]

with open('data/hg19_reference/chr'+chrom+'_oneline.txt') as infile:
    refseq=infile.read()

with open('data/hg19_chimp_align/human_chimp_diffs_chr'+chrom+'.txt') as infile:
    infile.next()
    anc_lines = [line.split() for line in infile if 'SNP' in line]

chimp_alleles = {}
for position, _, _, chimp_allele in anc_lines:
    chimp_alleles[int(position)] = chimp_allele

bases = 'ACGT'
mutations=[]
for trimer in product(bases, bases, bases):
    for base in bases:
        if trimer[1] != base:
            mutations.append((''.join(trimer), base))

print 'opening file'
infile=gzip.open('data/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
print 'file open'

for line in infile:
    if line.startswith('#CHROM'):
        break
print 'fast forwarded through file'

sample_ids = line.split()[9:]
n_lineages = 2 * len(sample_ids)

mutation_counts = {
    (mutation, haplotype_index): 0
    for haplotype_index in range(n_lineages)
    for mutation in mutations
}

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

    assert chrom == chromosome_number
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

anc_ind=0
for counter, line in enumerate(infile):
    try:
        (
            reference_allele,
            alternate_allele,
            position,
            derived_count,
            alleles
        ) = parse_line(line)
    except BadDataQualityError:
        continue

    context = refseq[position - 2 : position + 1]

    if 'N' not in context:
        if chimp_alleles.get(position) == alternate_allele:
            reverse=True
            derived_allele='0'
            this_mut=(context[0]+alternate_allele+context[2],reference_allele)
        else:
            reverse=False
            derived_allele='1'
            this_mut=(context,alternate_allele)

        if derived_count>1 and n_lineages-derived_count>1:
            if reverse:
                derived_count=n_lineages-derived_count
            i=9
            der_observed=0

            for i, allele in enumerate(alleles):
                if allele == derived_allele:
                    mutation_counts[(this_mut, i)] += 1
                    der_observed += 1
            assert der_observed == derived_count

    if counter > 100:
        break

def write_output(mutation_counts):
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

    with open('derived_each_lineage_chr'+chrom+'_nosingle.txt','w') as outfile:
        outfile.write(output)

write_output(mutation_counts)

print 'finished chrom ',chrom

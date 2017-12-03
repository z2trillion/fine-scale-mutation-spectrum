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
    anc_lines=infile.readlines()

bases = 'ACGT'
mutations=[]
for trimer in product(bases, bases, bases):
    for base in bases:
        if trimer[1] != base:
            mutations.append((''.join(trimer), base))

print 'opening file'
infile=gzip.open('data/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
print 'file open'

line=infile.readline()
while not line.startswith('#CHROM'):
    line=infile.readline()

print 'fast forwarded through file'
s=line.strip('\n').split('\t')
n_lineages=2*(len(s)-9)

mutation_counts = {
    (mutation, haplotype_index): 0
    for haplotype_index in range(n_lineages)
    for mutation in mutations
}


output='Mut_type'
for i in range(9,len(s)):
    output+=' '+sample_id_to_population[s[i]]
output+='\n'

popul=dict({})


indices = {group: [] for group in groups}

for i in range(9,len(s)):
    popul[i]=sample_id_to_population[s[i]]
    indices[population_to_group[popul[i]]].append(i)

anc_lines.pop(0)
anc_ind=0
while anc_ind<len(anc_lines):
    s=anc_lines[anc_ind].strip('\n').split(' ')
    if s[1]=='SNP':
        anc_lines[anc_ind]=deepcopy(s)
        anc_ind+=1
    else:
        anc_lines.pop(anc_ind)
anc_ind=0

for counter, line in enumerate(infile):
    s=line.strip('\n').split('\t')
    pos=int(s[1])
    context=refseq[pos-2:pos+1]
    if len(s[3]+s[4]+s[3]+s[4])==4 and s[6]=='PASS' and s[3] in 'ACGT' and s[3] in 'ACGT' and s[4] in 'ACGT' and s[4] in 'ACGT' and not 'N' in context:
        while anc_ind<len(anc_lines)-1 and int(anc_lines[anc_ind][0])<pos:
            anc_ind+=1
        if int(anc_lines[anc_ind][0])==pos and s[4]==anc_lines[anc_ind][3]:
            reverse=True
            der_allele='0'
            this_mut=(context[0]+s[4]+context[2],s[3])
        else:
            reverse=False
            der_allele='1'
            this_mut=(context,s[4])
        s2=s[7].split(';')
        count_der=int(s2[0][3:])
        if count_der>1 and n_lineages-count_der>1:
            if reverse:
                count_der=n_lineages-count_der
            i=9
            der_observed=0
            while i<len(s) and der_observed<count_der:
                for j in [0,2]:
                    if s[i][j]==der_allele:
                        mutation_counts[(this_mut,2*(i-9)+j/2)]+=1
                        der_observed+=1
                i+=1

    if counter > 100:
        break

for mut in mutations:
    output+=mut[0]+'_'+mut[1]
    for i in range(n_lineages):
        output+=' '+str(mutation_counts[(mut,i)])
    output+='\n'

outfile=open('derived_each_lineage_chr'+chrom+'_nosingle.txt','w')
outfile.write(output)
outfile.close()

print 'finished chrom ',chrom

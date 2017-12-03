from copy import deepcopy
import sys
import gzip
from itertools import product

chrom=sys.argv[1]

infile=open('data/1000genomes_phase3_sample_IDs.txt')
lines=infile.readlines()
infile.close()

big_populs=['EAS','SAS','EUR','AMR','AFR']
populs=dict({})
populs['EAS']=['CHB','JPT','CHS','CDX','KHV','CHD']
populs['EUR']=['CEU','TSI','GBR','FIN','IBS']
populs['AFR']=['YRI','LWK','GWD','MSL','ESN','ACB','ASW']
populs['SAS']=['GIH','PJL','BEB','STU','ITU']
populs['AMR']=['CLM','MXL','PUR','PEL']
group=dict({})
allpops=[]
for bigpop in big_populs:
    for pop in populs[bigpop]:
        allpops.append(pop)
        group[pop]=bigpop

pop_thisID=dict({})
for line in lines:
    s=line.split('\t')
    if len(s)>3:
        pop_thisID[s[0]]=s[2]


with open('data/hg19_reference/chr'+chrom+'_oneline.txt') as infile:
    refseq=infile.read()

with open('data/hg19_chimp_align/human_chimp_diffs_chr'+chrom+'.txt') as infile:
    infile.next()
    anc_lines = [line.split() for line in infile if 'SNP' in line]

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
#    print line
    line=infile.readline()

print 'fast forwarded through file'
s=line.strip('\n').split('\t')
num_lineages=2*(len(s)-9)

output='Mut_type'
for i in range(9,len(s)):
    output+=' '+pop_thisID[s[i]]
output+='\n'

popul=dict({})
indices=dict({})

for pop in allpops:
    indices[pop]=[]

for i in range(9,len(s)):
    popul[i]=pop_thisID[s[i]]
    indices[popul[i]].append(i)

count=dict({})

anc_ind=0

output=dict({})
mut_count=dict({})
for pop in allpops:
    output[pop]='Mut'
    for i in range(1,2*len(indices[pop])+1):
        output[pop]+=' '+str(i)
        for mut in mutations:
            mut_count[(mut,pop,i)]=0
    output[pop]+='\n'

for line in infile:
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
        if min(count_der,num_lineages-count_der)>1:
            if reverse:
                count_der=num_lineages-count_der
            i=9
            der_observed=0
            for pop in allpops:
                count[pop]=0
            while i<len(s) and der_observed<count_der:
                for j in [0,2]:
                    if s[i][j]==der_allele:
                        count[popul[i]]+=1
                        der_observed+=1
                i+=1
            for pop in allpops:
                if count[pop]>0:
                    mut_count[(this_mut,pop,count[pop])]+=1

for pop in allpops:
    for mut in mutations:
        output[pop]+=mut[0]+'_'+mut[1]
        for i in range(1,2*len(indices[pop])+1):
            output[pop]+=' '+str(mut_count[(mut,pop,i)])
        output[pop]+='\n'
    outfile=open('finescale_mut_spectra/mut_type_v_allele_freq_'+pop+'_chr'+chrom+'_nosingle.txt','w')
    outfile.write(output[pop])
    outfile.close()

print 'finished chrom ',chrom

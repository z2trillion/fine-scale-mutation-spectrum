from copy import deepcopy
import sys
import gzip

chrom=sys.argv[1]

infile=open('1000genomes_phase3_sample_IDs.txt')
lines=infile.readlines()
infile.close()

big_populs=['EAS','SAS','EUR','AMR','AFR']
populs=dict({})
populs['EAS']=['CHB','JPT','CHS','CDX','KHV','CHD']
populs['EUR']=['CEU','TSI','GBR','FIN','IBS']
populs['AFR']=['YRI','LWK','GWD','MSL','ESN']
populs['SAS']=['GIH','PJL','BEB','STU','ITU']
populs['AMR']=['CLM','MXL','PUR','PEL','ACB','ASW']
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


infile=open('../hg19_reference/chr'+chrom+'_oneline.txt')
refseq=infile.read()
infile.close()

infile=open('../hg19_chimp_align/human_chimp_diffs_chr'+chrom+'.txt')
anc_lines=infile.readlines()
infile.close()


muts=[]
for b1 in 'ACGT':
    for b2 in 'ACGT':
        for b3 in 'ACGT':
            for b4 in 'ACGT':
                if not b2==b4:
                    muts.append((b1+b2+b3,b4))

print 'opening file'
infile=gzip.open('../phase3_1kg/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz')
print 'file open'

line=infile.readline()
while not line.startswith('#CHROM'):
#    print line
    line=infile.readline()

print 'fast forwarded through file'
s=line.strip('\n').split('\t')
num_lineages=2*(len(s)-9)

mut_count=dict({})
for hap_ind in range(num_lineages):
    for mut in muts:
        mut_count[(mut,hap_ind)]=0

output='Mut_type'
for i in range(9,len(s)):
    output+=' '+pop_thisID[s[i]]
output+='\n'

popul=dict({})

indices=dict({})
for p in ['AMR','EUR','EAS','SAS','AFR']:
    indices[p]=[]


for i in range(9,len(s)):
    popul[i]=pop_thisID[s[i]]
    indices[group[popul[i]]].append(i)

count=dict({})
AF=dict({})

groups=['EUR','EAS','SAS','AMR','AFR']
group_info_indices=[8,5,9,6,7]

indcount=dict({})
indcount['0|0']=0
for gt in ['0|1','1|0','0/1','1/0']:
    indcount[gt]=1
indcount['1|1']=2

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
        if count_der>1 and num_lineages-count_der>1:
            if reverse:
                count_der=num_lineages-count_der
            i=9
            der_observed=0
            while i<len(s) and der_observed<count_der:
                for j in [0,2]:
                    if s[i][j]==der_allele:
                        mut_count[(this_mut,2*(i-9)+j/2)]+=1
                        der_observed+=1
                i+=1

for mut in muts:
    output+=mut[0]+'_'+mut[1]
    for i in range(num_lineages):
        output+=' '+str(mut_count[(mut,i)])
    output+='\n'

outfile=open('derived_each_lineage_chr'+chrom+'_nosingle.txt','w')
outfile.write(output)
outfile.close()

print 'finished chrom ',chrom


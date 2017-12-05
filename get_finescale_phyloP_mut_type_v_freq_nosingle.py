import sys
import gzip
from mutations import mutations, bases
from labels import sample_id_to_population, populations
from common import reference_sequence, human_chimp_differences

def get_finescale(dataset, chrom):
    if dataset == 'phastcons':
        infile_path = 'data/phastConsElements100way.txt'
        outfile_path = 'finescale_mut_spectra/phyloP_conserved_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
    elif dataset == 'nestedrepeats':
        infile_path = 'data/nestedRepeats.txt'
        outfile_path = 'finescale_mut_spectra/inrepeats_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'
    else:
        raise ValueError

    infile=open(infile_path)
    lines=infile.readlines()
    infile.close()

    ind=0
    s=lines[ind].split('\t')
    while not s[1]=='chr'+chrom:
        ind+=1
        s=lines[ind].split('\t')

    conserved=[(int(s[2]),int(s[3]))]
    ind+=1
    s=lines[ind].split('\t')

    print ind, len(lines),s
    while ind<len(lines)-1 and s[1]== 'chr'+chrom:
        if int(s[2])==conserved[-1][-1]+1:
            new_tup=(conserved[-1][0],int(s[3]))
            conserved.pop()
            conserved.append(new_tup)
        else:
            conserved.append((int(s[2]),int(s[3])))
        ind+=1
        line=lines[ind]
        s=line.strip('\n').split('\t')

    print len(conserved), conserved[:10]

    infile=open('data/1000genomes_phase3_sample_IDs.txt')
    lines=infile.readlines()
    infile.close()

    refseq = reference_sequence(chrom)
    anc_lines = human_chimp_differences(chrom)

    mut_count = {}
    for pop in populations:
        for mut in mutations:
            for i in range(250):
                mut_count[(mut,pop,i)]=0

    print 'opening file'
    infile=gzip.open('data/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
    print 'file open'

    line=infile.readline()
    while not line.startswith('#CHROM'):
        line=infile.readline()

    s=line.strip('\n').split('\t')
    num_lineages=2*(len(s)-9)

    output=dict({})
    for bigpop in populations:
        output[bigpop]='Ref Alt '

    popul=dict({})

    indices=dict({})
    for p in populations:
        indices[p]=[]

    for i in range(9,len(s)):
        popul[i]=sample_id_to_population[s[i]]
        indices[popul[i]].append(i)

    count = {}
    AF=dict({})

    anc_ind=0

    conserved_ind=0

    for line_counter, line in enumerate(infile):
        s=line.strip('\n').split('\t')
        pos=int(s[1])
        context=refseq[pos-2:pos+1]
        while conserved_ind<len(conserved)-1 and pos>conserved[conserved_ind][1]:
            conserved_ind+=1
        if pos>= conserved[conserved_ind][0] and pos<=conserved[conserved_ind][1] and len(s[3]+s[4])==2 and s[6]=='PASS' and s[3] in 'ACGT' and s[4] in 'ACGT' and not 'N' in context:
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
                for pop in populations:
                    count[pop]=0
                while i<len(s) and der_observed<count_der:
                    for j in [0,2]:
                        if s[i][j]==der_allele:
                            count[popul[i]]+=1
                            der_observed+=1
                    i+=1
                for pop in populations:
                    if count[pop]>0:
                        mut_count[(this_mut,pop,count[pop])]+=1
        print line_counter
        if line_counter > 1e4:
            break

    for pop in populations:
        for mut in mutations:
            output[pop]+=mut[0]+'_'+mut[1]
            for i in range(1,2*len(indices[pop])+1):
                output[pop]+=' '+str(mut_count[(mut,pop,i)])
            output[pop]+='\n'
        outfile=open(outfile_path %  pop,'w')
        outfile.write(output[pop])
        outfile.close()

    print 'finished chrom ',chrom

if __name__ == '__main__':
    chrom=sys.argv[1]
    get_finescale('phastcons', chrom)

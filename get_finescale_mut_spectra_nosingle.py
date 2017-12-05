import sys

from mutations import mutations, bases
from labels import sample_id_to_population, populations
from common import reference_sequence, human_chimp_differences

from common import (
    write_output,
    get_column_indices,
    get_column_index_to_population,
    initialize_mut_count,
    open_infile,
)

def get_finescale(chrom):
    outfile_path = 'finescale_mut_spectra/mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'

    infile, line = open_infile(chrom)
    s=line.strip('\n').split('\t')
    num_lineages=2*(len(s)-9)

    indices = get_column_indices(s)
    popul = get_column_index_to_population(s)
    mut_count = initialize_mut_count(indices)

    output={}

    for pop in populations:
        output[pop]='Mut'
        for i in range(1,2*len(indices[pop])+1):
            output[pop]+=' '+str(i)
        output[pop]+='\n'

    anc_ind=0
    refseq = reference_sequence(chrom)
    anc_lines = human_chimp_differences(chrom)
    for line_counter, line in enumerate(infile):
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
                count = {population: 0 for population in populations}
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

    write_output(output, outfile_path, indices, mut_count)

    print 'finished chrom ',chrom

if __name__ == '__main__':
    chrom=sys.argv[1]
    get_finescale(chrom)

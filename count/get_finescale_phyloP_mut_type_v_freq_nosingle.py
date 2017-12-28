import sys

from mutations import mutations, bases
from labels import sample_id_to_population, populations
from common import (
    reference_sequence,
    get_human_chimp_differences,
    write_output,
    get_column_indices,
    get_column_index_to_population,
    initialize_mut_count,
    open_infile,
    get_conserved,
)

def get_finescale(outfile_path, chrom, output, conserved):
    infile, line = open_infile(chrom)
    column_labels = line.split()
    n_lineages = 2 * (len(column_labels) - 9)

    indices = get_column_indices(column_labels)
    column_index_to_population = get_column_index_to_population(column_labels)
    mut_count = initialize_mut_count(indices)

    conserved_ind=0
    refseq = reference_sequence(chrom)
    human_chimp_differences = get_human_chimp_differences(chrom)
    for line_counter, line in enumerate(infile):
        s=line.strip('\n').split('\t')
        pos=int(s[1])
        context = refseq[pos-2:pos+1]
        if 'N' in context:
            continue
        while conserved_ind<len(conserved)-1 and pos>conserved[conserved_ind][1]:
            conserved_ind+=1
        if pos>= conserved[conserved_ind][0] and pos<=conserved[conserved_ind][1] and len(s[3]+s[4])==2 and s[6]=='PASS' and s[3] in 'ACGT' and s[4] in 'ACGT':
            if human_chimp_differences.get(pos) == s[4]:
                reverse=True
                der_allele='0'
                this_mut=(context[0]+s[4]+context[2],s[3])
            else:
                reverse=False
                der_allele='1'
                this_mut=(context,s[4])
            s2=s[7].split(';')
            count_der=int(s2[0][3:])
            if min(count_der,n_lineages-count_der)>1:
                if reverse:
                    count_der=n_lineages-count_der
                i=9
                der_observed=0
                count = {population: 0 for population in populations}
                while i<len(s) and der_observed<count_der:
                    for j in [0,2]:
                        if s[i][j]==der_allele:
                            count[column_index_to_population[i]]+=1
                            der_observed+=1
                    i+=1
                for pop in populations:
                    if count[pop]>0:
                        mut_count[(this_mut,pop,count[pop])]+=1
        # print line_counter
        if line_counter > 1e4:
            break

    write_output(output, outfile_path, indices, mut_count)

    print 'finished chrom ',chrom

if __name__ == '__main__':
    chrom=sys.argv[1]
    output = {population: 'Ref Alt \n' for population in populations}
    outfile_path = '../finescale_mut_spectra/phyloP_conserved_mut_type_v_allele_freq_%s_chr'+chrom+'_nosingle.txt'

    conserved = get_conserved('../data/bed_files/phastConsElements100way.txt', chrom)
    get_finescale(outfile_path, chrom, output, conserved)

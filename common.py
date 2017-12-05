from mutations import mutations, bases
from labels import sample_id_to_population, populations
import gzip

def reference_sequence(chromosome_number):
    with open('data/hg19_reference/chr'+chromosome_number+'_oneline.txt') as infile:
        return infile.read()

def human_chimp_differences(chromosome_number):
    with open('data/hg19_chimp_align/human_chimp_diffs_chr'+chromosome_number+'.txt') as infile:
        infile.next()

        return [line.split() for line in infile if 'SNP' in line]

def get_column_indices(column_labels):
    population_to_column_indices = {population: [] for population in populations}

    sample_ids = column_labels[9:]
    for i, sample_id in enumerate(sample_ids):
        population_to_column_indices[sample_id_to_population[sample_id]].append(i + 9)

    return population_to_column_indices

def get_column_index_to_population(column_labels):
    sample_ids = column_labels[9:]
    return {
        i + 9 : sample_id_to_population[sample_id]
        for i, sample_id in enumerate(sample_ids)
    }

def initialize_mut_count(population_to_column_indices):
    mut_count = {}
    for pop in populations:
        for mut in mutations:
            for i in range(1,2*len(population_to_column_indices[pop])+1):
                mut_count[(mut,pop,i)]=0
    return mut_count


def write_output(output, outfile_path, indices, mut_count):
    for pop in populations:
        for mut in mutations:
            output[pop]+=mut[0]+'_'+mut[1]
            for i in range(1,2*len(indices[pop])+1):
                output[pop]+=' '+str(mut_count[(mut,pop,i)])
            output[pop]+='\n'
        outfile=open(outfile_path %  pop,'w')
        outfile.write(output[pop])
        outfile.close()

def open_infile(chrom):
    print 'opening file'
    infile=gzip.open('data/ALL.chr'+chrom+'.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz')
    print 'file open'

    line=infile.readline()
    while not line.startswith('#CHROM'):
        line=infile.readline()

    print 'fast forwarded through file'

    return infile, line

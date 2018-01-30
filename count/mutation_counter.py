from mutations import mutations, bases
from labels import sample_id_to_population, populations
from common import (
    reference_sequence,
    get_human_chimp_differences,
    write_output,
    get_column_indices,
    get_column_index_to_population,
    initialize_mut_count,
)

class MutationCounter:
    def __init__(self, chrom, included_regions):
        self.chrom = chrom
        self.included_regions = included_regions

    def configure(self, line):
        column_labels = line.split()
        self.n_lineages = 2 * (len(column_labels) - 9)
        self.column_index_to_population = get_column_index_to_population(column_labels)

        self.refseq = reference_sequence(self.chrom)
        self.human_chimp_differences = get_human_chimp_differences(self.chrom)

        for included_region in self.included_regions:
            included_region.configure(column_labels)

    def process_line(self, line):
        s=line.strip('\n').split('\t')
        pos=int(s[1])
        context = self.refseq[pos-2:pos+1]
        if 'N' in context:
            return

        for included_region in self.included_regions:
            included_region.update_position(pos)

        if len(s[3]+s[4])==2 and s[6]=='PASS' and s[3] in 'ACGT' and s[4] in 'ACGT':
            if self.human_chimp_differences.get(pos) == s[4]:
                reverse=True
                der_allele='0'
                this_mut=(context[0]+s[4]+context[2],s[3])
            else:
                reverse=False
                der_allele='1'
                this_mut=(context,s[4])
            s2=s[7].split(';')
            count_der=int(s2[0][3:])
            if min(count_der,self.n_lineages-count_der)>1:
                if reverse:
                    count_der=self.n_lineages-count_der
                i=9
                der_observed=0
                count = {population: 0 for population in populations}
                while i<len(s) and der_observed<count_der:
                    for j in [0,2]:
                        if s[i][j]==der_allele:
                            count[self.column_index_to_population[i]]+=1
                            der_observed+=1
                    i+=1
                for pop in populations:
                    if count[pop]>0:
                        for included_regions in self.included_regions:
                            included_region.update_counts(pos, this_mut, pop, count[pop])

    def write_output(self):
        for included_region in self.included_regions:
            included_region.write_output()


class IncludedRegion:
    def __init__(self, chrom, output, outfile_path, conserved):
        self.chrom = chrom
        self.output = output
        self.outfile_path = outfile_path
        self.conserved = conserved

        self.conserved_ind = 0

    def configure(self, column_labels):
        self.indices = get_column_indices(column_labels)
        self.column_index_to_population = get_column_index_to_population(column_labels)
        self.mut_count = initialize_mut_count(self.indices)

    def update_position(self, pos):
        while self.conserved_ind<len(self.conserved)-1 and pos>self.conserved[self.conserved_ind][1]:
            self.conserved_ind+=1

    def update_counts(self, position, mutation, population, count):
        if (
            position >= self.conserved[self.conserved_ind][0] and
            position <= self.conserved[self.conserved_ind][1]
        ):
            self.mut_count[(mutation, population, count)] += 1

    def write_output(self):
        write_output(self.output, self.outfile_path, self.indices, self.mut_count)

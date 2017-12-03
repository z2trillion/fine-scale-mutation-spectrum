
def reference_sequence(chromosome_number):
    with open('data/hg19_reference/chr'+chromosome_number+'_oneline.txt') as infile:
        return infile.read()

def human_chimp_differences(chromosome_number):
    with open('data/hg19_chimp_align/human_chimp_diffs_chr'+chromosome_number+'.txt') as infile:
        infile.next()

        return [line.split() for line in infile if 'SNP' in line]

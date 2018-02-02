import os
# import gzip


def gen_oneline(chrom):
    # infile=gzip.open('../data/hg19_reference/chr%i.fa.gz' % chrom)
    infile = open('../data/hg19_reference/chr%i.fa' % chrom)
    lines = infile.readlines()
    infile.close()

    output = ''
    for line in lines[1:]:
        output += line.upper().strip('\n')

    outfile = open('../data/hg19_reference/chr%i_oneline.txt' % chrom, 'w')
    outfile.write(output)
    outfile.close()


if __name__ == '__main__':
    for chrom in range(1, 23) + ['X']:
        if os.path.isfile('../data/hg19_reference/chr%i_oneline.txt' % chrom):
            print 'existing oneline file found for chromosome %i' % chrom
        else:
            gen_oneline(chrom)
            print 'generated oneline file for chromosome %i' % chrom

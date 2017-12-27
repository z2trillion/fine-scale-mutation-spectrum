import sys

chrom=sys.argv[1]

infile=open('chr'+chrom+'.fa')
lines=infile.readlines()
infile.close()

output=''
for line in lines[1:]:
    output+=line.upper().strip('\n')

outfile=open('chr'+chrom+'_oneline.txt','w')
outfile.write(output)
outfile.close()

import os
import gzip

def extract_diffs(chrom):
    infile=gzip.open('../data/hg19_chimp_align/chr%i.hg19.panTro4.net.axt.gz' % chrom)
    lines=infile.readlines()
    infile.close()

    line_ind=0
    while not lines[line_ind][0]=='0':
        line_ind+=1

    output='Pos SNP/Indel Human Chimp\n'

    last_bin_start=0
    last_end=0
    while line_ind<len(lines):
        s=lines[line_ind].split(' ')
    #    print s
        start, end = int(s[2]), int(s[3])
        if start>last_end+1:
            output+=str(last_end)+' Indel '+str(start-last_end-1)+'\n'
        last_end=end
        if start>last_bin_start+10**6:
            # print start
            last_bin_start=start
        human=lines[line_ind+1].upper().strip('\n')
        chimp=lines[line_ind+2].upper().strip('\n')
        char_ind=0
        while char_ind<len(human):
            found_indel=False
            h_allele=''
            c_allele=''
    #        print human[char_ind], chimp[char_ind]
            if human[char_ind]==chimp[char_ind] or human[char_ind]=='N' or chimp[char_ind]=='N':
                char_ind+=1
            elif human[char_ind] in 'ACGT' and chimp[char_ind] in 'ACGT':
                found_match=False
                h_allele=human[char_ind]
                c_allele=chimp[char_ind]
                char_ind+=1
                while not found_match and not found_indel and char_ind<len(human):
                    if human[char_ind]==chimp[char_ind]:
                        found_match=True
                    elif human[char_ind]=='-' or chimp[char_ind]=='-':
                        found_indel=True
                    else:
                        h_allele+=human[char_ind]
                        c_allele+=chimp[char_ind]
                        char_ind+=1
                if not found_indel:
                    for i in range(len(h_allele)):
                        j=char_ind-len(h_allele)+i
                        output+=str(start+j)+' SNP '+human[j]+' '+chimp[j]+'\n'
            if char_ind<len(human) and (human[char_ind]=='-' or chimp[char_ind]=='-'):
                output+=str(start+char_ind)+' Indel '
                len_indel=0
                while char_ind<len(human) and not human[char_ind]==chimp[char_ind]:
                    if human[char_ind]=='-':
                        start-=1
    #                h_allele+=human[char_ind]
    #                c_allele+=chimp[char_ind]
                    len_indel+=1
                    char_ind+=1
                output+=str(len_indel)+'\n'
    #            output+=h_allele+' '+c_allele+'\n'
    #        if not h_allele=='':
    #            print h_allele, c_allele, end-(start+char_ind)
        line_ind+=4

    output+=str(last_end+1)+' Indel 100000000\n'

    outfile=open('../data/hg19_chimp_align/human_chimp_diffs_chr%i.txt' % chrom,'w')
    outfile.write(output)
    outfile.close()

if __name__ == '__main__':
    for chrom in range(1, 23) + ['X']:
        if os.path.isfile('../data/hg19_chimp_align/human_chimp_diffs_chr%i.txt' % chrom):
            print 'existing human_chimp_diffs file found for chromosome %i' % chrom
        else:
            extract_diffs(chrom)
            print 'extracted diffs for chromosome %i' % chrom

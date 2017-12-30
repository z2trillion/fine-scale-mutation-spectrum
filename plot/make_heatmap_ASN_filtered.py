import sys
import matplotlib
matplotlib.use('Agg')  # This prevents the plotting engine from starting up.
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import chi2_contingency
import math

comp=dict({})
comp['A'],comp['C'],comp['G'],comp['T']='T','G','C','A'
ypos, ylabel=[],[]

inv_mut_index=dict({})
mut_index=dict({})
row, col = 0,0

for (b2,d) in [('A','T'),('A','C'),('A','G'),('C','T'),('C','G'),('C','A')]:
    for b1 in 'ACGT':
        col=0
        ypos.append(row+0.5)
        if b1=='T' and b2=='C' and d=='A':
            ylabel.append('5\'-'+b1)
        elif b1=='C':
            ylabel.append(b2+r'$\to$'+d+r'  '+b1)
        else:
            ylabel.append(b1)
        for b3 in 'ACGT':
            mut_index[(b1+b2+b3,d)]=(row,col)
            inv_mut_index[(row,col)]=b1+b2+b3+'_'+d
            mut_index[(comp[b3]+comp[b2]+comp[b1],comp[d])]=(row,col)
            col+=1
        row+=1

def frequency_breakdown(pop,start_chr):
    count_array=np.zeros((row,col))
    for chrom in range(start_chr,5):
        infile=open('../finescale_mut_spectra/mut_type_v_allele_freq_'+pop+'_chr'+str(chrom)+'_nosingle.txt')
        lines=infile.readlines()
        infile.close()

        s=lines[1].strip('\n').split(' ')
        start_ind=2
        end_ind=len(s)-2
        while 1.0*end_ind/(len(s)-2)>0.98:
            end_ind-=1
        for line in lines[1:]:
            s=line.strip('\n').split(' ')
            for i in range(start_ind,end_ind):
                count_array[mut_index[(s[0][:3],s[0][4])]]+=int(s[i])
    return count_array

def frequency_breakdown_phyloP(pop,start_chr):
    count_array=np.zeros((row,col))
    for chrom in range(start_chr,5):
        infile=open('../finescale_mut_spectra/phyloP_conserved_mut_type_v_allele_freq_'+pop+'_chr'+str(chrom)+'_nosingle.txt')
        lines=infile.readlines()
        infile.close()

        s=lines[0].strip('\n').split(' ')
        start_ind=2
        end_ind=len(s)-1
        while 1.0*end_ind/(len(s)-2)>0.98:
            end_ind-=1
        for line in lines[1:]:
            s=line.strip('\n').split(' ')
            for i in range(start_ind,end_ind):
                count_array[mut_index[(s[0][:3],s[0][4])]]+=int(s[i])
    return count_array

def frequency_breakdown_repeats(pop,start_chr):
    count_array=np.zeros((row,col))
    for chrom in range(start_chr,5):
        infile=open('../finescale_mut_spectra/inrepeats_mut_type_v_allele_freq_'+pop+'_chr'+str(chrom)+'_nosingle.txt')
        lines=infile.readlines()
        infile.close()

        s=lines[0].strip('\n').split(' ')
        start_ind=2
        end_ind=len(s)-1
        while 1.0*end_ind/(len(s)-2)>0.98:
            end_ind-=1
        for line in lines[1:]:
            s=line.strip('\n').split(' ')
            for i in range(start_ind,end_ind):
                count_array[mut_index[(s[0][:3],s[0][4])]]+=int(s[i])
    return count_array

pop_counts=dict({})
num_variants=dict({})
pops=['JPT','CHB','CHS','CDX','KHV']
for pop in pops:
    pop_counts[pop]=frequency_breakdown(pop,1)-frequency_breakdown_repeats(pop,1)-frequency_breakdown_phyloP(pop,1)
    num_variants[pop]=pop_counts[pop].sum()

name=dict({})
for p in zip(pops,['Japanese','Northern Han','Southern Han','Dai','Kinh']):
    name[p[0]]=p[1]

for pop_ind1 in range(5):
    refpop=pops[pop_ind1]
    subplot_ind=1
    for pop_ind2 in range(pop_ind1)+range(pop_ind1+1,5):
        ratio_list=[]
        pop=pops[pop_ind2]
        print refpop, pop
        ratio_grid=np.zeros((row,col))
        sig_x,sig_y=[],[]
        for i in range(row):
            for j in range(col):
                chi2_results=chi2_contingency(np.array([[pop_counts[pop][i][j],num_variants[pop]],[pop_counts[refpop][i][j],num_variants[refpop]]]))
                this_pval=chi2_results[1]
                ratio_grid[i][j]=pop_counts[pop][i][j]*num_variants[refpop]/(num_variants[pop]*pop_counts[refpop][i][j])
                if this_pval<0.00001:
                    sig_x.append(j+0.5)
                    sig_y.append(i+0.5)
        plt.subplot(1,5,subplot_ind)
        subplot_ind+=1
        plt.title(name[pop]+' v\n'+name[refpop],fontsize=10)

        if subplot_ind==2:
            plt.xticks((0.5,1.5,2.5,3.5),('3\'-A','C','G','T'))
            plt.yticks(tuple(ypos),tuple(ylabel))
        else:
            plt.xticks((0.5,1.5,2.5,3.5),('A','C','G','T'))
            plt.yticks(())
        for k in range(1,6):
            plt.axhline(y=k*4,color='black')
        plt.pcolor(ratio_grid,vmin=0.84,vmax=1.16,cmap='seismic')
        if subplot_ind==5:
            plt.colorbar()
        plt.scatter(sig_x,sig_y,marker='.')

    fig=plt.gcf()
    plt.savefig('ASN_heatmap_v_'+refpop+'_nosingle.pdf',format='pdf')
    plt.clf()

import sys
import argparse
import numpy as np
from scipy.stats import chi2_contingency
from itertools import product

import matplotlib
matplotlib.use('Agg')  # This prevents the plotting engine from starting up.
import matplotlib.pyplot as plt  # noqa E402

comp = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
}
ypos, ylabel = [], []

inv_mut_index = {}
mut_index = {}
row, col = 0, 0

for b2, d in [('A', 'T'), ('A', 'C'), ('A', 'G'),
              ('C', 'T'), ('C', 'G'), ('C', 'A')]:
    for b1 in 'ACGT':
        col = 0
        ypos.append(row+0.5)
        if b1 == 'T' and b2 == 'C' and d == 'A':
            ylabel.append('5\'-'+b1)
        elif b1 == 'C':
            ylabel.append(b2+r'$\to$'+d+r'  '+b1)
        else:
            ylabel.append(b1)
        for b3 in 'ACGT':
            mut_index[(b1+b2+b3, d)] = (row, col)
            inv_mut_index[(row, col)] = b1+b2+b3+'_'+d
            mut_index[(comp[b3]+comp[b2]+comp[b1], comp[d])] = (row, col)
            col += 1
        row += 1


def frequency_breakdown(path, chromosomes, frequency_range):
    count_array = np.zeros((row, col))
    for chrom in chromosomes:
        infile = open(path % str(chrom))
        lines = infile.readlines()
        infile.close()

        s = lines[1].strip('\n').split(' ')
        # TODO(mason) do these better.
        start_index = 2
        while float(start_index - 2) / (len(s) - 2) < frequency_range[0]:
            start_index += 1

        end_index = len(s) - 1
        while float(end_index) / (len(s)-2) > frequency_range[1]:
            end_index -= 1

        for line in lines[1:]:
            s = line.strip('\n').split(' ')
            for i in range(start_index, end_index):
                count_array[mut_index[(s[0][:3], s[0][4])]] += int(s[i])
    return count_array


def heatmap(chromosomes, population_pair, frequency_range, exclude, p_value):
    pop_counts = {}
    num_variants = {}

    for pop in population_pair:
        path = ('../finescale_mut_spectra/mut_type_v_allele_freq_' +
                pop + '_chr%s_nosingle.txt')
        pop_counts[pop] = frequency_breakdown(path, chromosomes,
                                              frequency_range)
        if exclude:
            repeats_path = ('../finescale_mut_spectra/'
                            'inrepeats_mut_type_v_allele_freq_' +
                            pop + '_chr%s_nosingle.txt')
            pop_counts[pop] -= frequency_breakdown(repeats_path, chromosomes,
                                                   frequency_range)

            conserved_path = ('../finescale_mut_spectra/'
                              'phyloP_conserved_mut_type_v_allele_freq_' +
                              pop + '_chr%s_nosingle.txt')
            pop_counts[pop] -= frequency_breakdown(conserved_path, chromosomes,
                                                   frequency_range)
        num_variants[pop] = pop_counts[pop].sum()

    refpop, pop = population_pair

    ratio_grid = np.zeros((row, col))
    sig_x, sig_y = [], []
    for i in range(row):
        for j in range(col):
            _, this_pval, _, _ = chi2_contingency(
                np.array([
                    [pop_counts[pop][i][j], num_variants[pop]],
                    [pop_counts[refpop][i][j], num_variants[refpop]]
                ])
            )
            ratio_grid[i][j] = (pop_counts[pop][i][j] * num_variants[refpop] /
                                (num_variants[pop] * pop_counts[refpop][i][j]))
            if this_pval < p_value:
                sig_x.append(j+0.5)
                sig_y.append(i+0.5)

    return ratio_grid, (sig_x, sig_y)


def make_titles(chromosome_groups, population_pairs, frequency_range, exclude,
                p_value):
    title = ''

    distinct_population_pairs = len(set(population_pairs)) != 1
    if not distinct_population_pairs:
        title += '/'.join(population_pairs[0])

    distinct_chromosome_groups = (
        len(set(tuple(cg) for cg in chromosome_groups)) != 1
    )
    if not distinct_chromosome_groups:
        if set(chromosome_groups[0]) == set(str(c) for c in range(1, 23)):
            title += ' Autosomes'
        elif len(chromosome_groups[0]) == 1:
            title += ' Chromosome %s' % str(chromosome_groups[0][0])
        else:
            title += ' Chromosomes %s' % ' '.join(
                str(c) for c in chromosome_groups[0]
            )
    # if exclude:
    #     title += '\nExcluding Repeats and Conserved Regions'
    if frequency_range != [0, 1]:
        title += ' %.2f <= f <= %.2f' % tuple(frequency_range)
    title += ' p < %1.1e' % p_value

    column_titles = []
    for chromosomes, population_pair in product(chromosome_groups,
                                                population_pairs):
        column_title = ''
        if distinct_population_pairs:
            column_title += '/'.join(population_pair)
        if distinct_chromosome_groups:
            column_title += ' '
            column_title += ' '.join(str(c) for c in chromosomes)
        column_titles.append(column_title)

    return title, column_titles


def make_plot(ratio_grids, significant_indices, plot_title, column_titles):
    combined_grids = np.hstack(ratio_grids)
    r = max(1-np.min(combined_grids), np.max(combined_grids)-1)
    plt.pcolor(np.hstack(ratio_grids), vmin=1-r, vmax=1+r, cmap='seismic')
    plt.colorbar()

    xtick_values = len(ratio_grids) * list('ACGT')
    xtick_values[0] = "3'-A"
    plt.xticks(np.arange(.5, 4 * len(ratio_grids)), xtick_values, rotation=90)
    plt.yticks(ypos, ylabel)

    for k in range(len(ratio_grids)):
        plt.axvline(x=4 * k, color='black')

    for k in range(0, len(ratio_grids[0]), 4):
        plt.axhline(y=k, color='black', linestyle='dashed')

    top_x_axis = plt.twiny()  # Create a twin Axes sharing the yaxis, lol
    top_x_axis.set_xticks(np.arange(2, 4 * len(ratio_grids), 4))
    top_x_axis.set_xticklabels(column_titles)

    combined_x_positions = []
    combined_y_positions = []
    for i, (x_positions, y_positions) in enumerate(significant_indices):
        for x, y in zip(x_positions, y_positions):
            combined_x_positions.append(x + 4 * i)
            combined_y_positions.append(y)
    plt.scatter(combined_x_positions, combined_y_positions, marker='.',
                color='white')

    plt.xlim(0, len(combined_grids[0]))
    plt.ylim(0, len(combined_grids))

    plt.title(plot_title + '\n')

    plt.gcf()
    plt.savefig('heatmap.pdf', format='pdf')
    plt.clf()


valid_populations = [
    'CHB', 'JPT', 'CHS', 'CDX', 'KHV', 'CHD',
    'CEU', 'TSI', 'GBR', 'FIN', 'IBS',
    'YRI', 'LWK', 'GWD', 'MSL', 'ESN', 'ACB', 'ASW',
    'GIH', 'PJL', 'BEB', 'STU', 'ITU',
    'CLM', 'MXL', 'PUR', 'PEL',
]


def Population(population):
    if population not in valid_populations:
        raise ValueError
    return population


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=[str(c) for c in range(1, 23)])
    parser.add_argument('-i', '--individually', action='store_true')
    parser.add_argument('-p', '--population-pairs', type=Population,
                        nargs='+')
    parser.add_argument('-f', '--frequency-range', type=float, nargs=2,
                        default=[0, 1])
    parser.add_argument('-e', '--exclude', action='store_true')
    parser.add_argument('--p-value', type=float, default=1e-5)

    args = parser.parse_args(sys.argv[1:])

    chromosomes = args.chromosomes
    for chromosome in chromosomes:
        assert chromosome == 'X' or int(chromosome) in range(1, 23)
    if args.individually:
        chromosome_groups = [[chromosome] for chromosome in chromosomes]
    else:
        chromosome_groups = [chromosomes]

    population_pairs = zip(args.population_pairs[::2],
                           args.population_pairs[1::2])
    frequency_range = args.frequency_range
    exclude = args.exclude
    p_value = args.p_value

    heatmaps = [
        heatmap(
            c, population_pair, frequency_range, exclude, p_value
        ) for c, population_pair in product(chromosome_groups,
                                            population_pairs)
    ]

    ratio_grids, significant_indices = zip(*heatmaps)

    plot_title, column_titles = make_titles(
        chromosome_groups, population_pairs, frequency_range, exclude, p_value
    )

    make_plot(ratio_grids, significant_indices, plot_title, column_titles)

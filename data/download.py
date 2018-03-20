import argparse
import errno
import os
import sys
import urllib


def get_file_name(chromosome):
    version_number = 1 if chromosome == 'X' else 5
    return ('ALL.chr%s.phase3_shapeit2_mvncall_integrated_'
            'v%ia.20130502.genotypes.vcf.gz' % (chromosome, version_number))


def download_vcf(chromosome):
    file_name = get_file_name(chromosome)
    file_path = os.path.join('vcfs', file_name)
    if not os.path.isfile(file_path):
        label = 'Downloading %s, %i%% complete'

        def progress_bar(block_number, block_size, file_size):
            if block_number % 100 == 0:
                fraction = 100 * block_number * block_size / file_size
                sys.stdout.write('\r' + label % (file_name, fraction))
                sys.stdout.flush()

        vcf_url = ('http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/'
                   'phase3/' + file_name)
        urllib.URLopener().retrieve(vcf_url, file_path, progress_bar)
        print


def download_reference_sequence(chromosome):
    file_name = 'chr%s.fa.gz' % chromosome
    file_path = os.path.join('hg19_reference', file_name)
    if not os.path.isfile(file_path):
        url = ('http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes/' +
               file_name)
        urllib.URLopener().retrieve(url, file_path)


def download_human_chimp_alignment(chromosome):
    file_name = 'chr%s.hg19.panTro4.net.axt.gz' % chromosome
    file_path = os.path.join('hg19_chimp_align', file_name)
    if not os.path.isfile(file_path):
        url = (
            'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/vsPanTro4/axtNet/'
            + file_name
        )
        urllib.URLopener().retrieve(url, file_path)


def download_sample_ids_and_bed_files():
    file_path_download_url_map = {
        '1000genomes_phase3_sample_IDs.txt': (
            'https://raw.githubusercontent.com/LukeAndersonTrocme/MutSpect/'
            'master/1000genomes_phase3_sample_IDs.txt'
        ),
        'bed_files/phastConsElements100way.txt.gz': (
            'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/'
            'phastConsElements100way.txt.gz'
        ),
        'bed_files/nestedRepeats.txt.gz': (
            'http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/'
            'nestedRepeats.txt.gz'
        )
    }

    for file_path, download_url in file_path_download_url_map.items():
        if not os.path.isfile(file_path):
            urllib.URLopener().retrieve(download_url, file_path)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--chromosomes', type=str, nargs='+',
                        default=['X'] + map(str, range(1, 23)))
    chromosomes = parser.parse_args(sys.argv[1:]).chromosomes

    for subdirectory in ['vcfs', 'hg19_reference', 'hg19_chimp_align',
                         'bed_files']:
        try:
            os.mkdir(subdirectory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    download_sample_ids_and_bed_files()

    for chromosome in chromosomes:
        download_reference_sequence(chromosome)
        download_vcf(chromosome)
        download_human_chimp_alignment(chromosome)

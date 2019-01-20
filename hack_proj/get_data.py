from utils import *
import sys

THRESHOLD = 500
SEQ_LEN = 1001


def get_locations(genome_tag_range):
    """
    :return: Two lists of lists with the location of cpg island and not
    """
    # d1 = genome_tag_range.get_untagged_ranges('chr1', distance=1000, range_len=2000)
    # d1 = genome_tag_range.get_random_untagged_ranges('chr1', ranges_num=100, range_len=10000)
    d1 = genome_tag_range.get_tagged_ranges('chr1')
    loc_cpg_island = list()
    for start, end in zip(d1[RANGE_START], d1[RANGE_END]):
        # print('%d: %d -- %d' % (count, start, end))
        loc_cpg_island.append([start, end])
    loc_not_cpg_island = list()
    d1 = genome_tag_range.get_untagged_ranges('chr1', THRESHOLD, SEQ_LEN)
    for start, end in zip(d1[RANGE_START], d1[RANGE_END]):
        # print('%d: %d -- %d' % (count, start, end))
        loc_not_cpg_island.append([start, end])

    return loc_cpg_island, loc_not_cpg_island


def get_chrom(chrom_name):
    """
    :param chrom: the wanted chromosome
    :return: string represents the chromosome
    """
    with open(chrom_name, 'r') as f:
        chrom_string = f.read()
    return chrom_string


def set_data(chrom_name, chrom_string, loc_cpg_island, loc_not_cpg_island, genome_tag_range):
    """
    :param chrom: the wanted chromosome as string
    :param loc_cpg_island: list of lists with the start point and end point of all the cpgisland
    :param loc_not_cpg_island: list of lists with the start point and end point of all the not cpg island
    :return: two lists: a list of the seqs and a list of their labels (accordingly)
    """
    seqs_list = list()
    labels_list = list()
    nucleotids_labels = list()
    for not_cpg in loc_not_cpg_island:
        seqs_list.append(chrom_string[not_cpg[0]:not_cpg[1]].upper())
        labels_list.append(0)
        start, length = not_cpg[0], not_cpg[1] - not_cpg[0]
        nucleotids_labels.append(genome_tag_range.get_tag_for_range(chrom_name, start, length))
    for cpg in loc_cpg_island:
        if cpg[1] - cpg[0] < THRESHOLD:
            continue
        seqs_list.append(chrom_string[cpg[0]:cpg[0] + SEQ_LEN].upper())
        labels_list.append(1)
        start, length = cpg[0], cpg[1] - cpg[0]
        nucleotids_labels.append(genome_tag_range.get_tag_for_range(chrom_name, start, length))

    return seqs_list, labels_list, nucleotids_labels


if __name__ == '__main__':
    """ sys.argv[1] = the wanted chromosome's name 
        sys.argv[2] = the location of the CGI.hg19.bed file
        sys.argv[3] = the location of the hg19.chrom.sizes file
    """
    chrom_name = sys.argv[1]
    genome_tag_range = GenomeTagRange(sys.argv[2], sys.argv[3])
    loc_cpg_island, loc_not_cpg_island = get_locations(genome_tag_range)
    chrom_string = get_chrom(chrom_name)
    seqs_list, labels_list, nucleotids_labels = set_data(chrom_name, chrom_string, loc_cpg_island, loc_not_cpg_island, genome_tag_range)


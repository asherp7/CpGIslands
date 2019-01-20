from utils import *
import sys

THRESHOLD = 500
SEQ_LEN = 1001


def get_locations():
    """
    :return: Two lists of lists with the location of cpg island and not
    """
    genome_tag_range = GenomeTagRange(
        r'C:\Users\User\Documents\Thirdyear\Tommy_course\CpGIslands\hack_proj\data\CGI.hg19.bed',
        r'C:\Users\User\Documents\Thirdyear\Tommy_course\CpGIslands\hack_proj\data\hg19.chrom.sizes')
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


def get_chrom(chrom):
    """
    :param chrom: the wanted chromosome
    :return: string represents the chromosome
    """
    with open(chrom, 'r') as f:
        chrom = f.read()
    return chrom


def set_data(chrom, loc_cpg_island, loc_not_cpg_island):
    """
    :param chrom: the wanted chromosome
    :param loc_cpg_island: list of lists with the start point and end point of all the cpgisland
    :param loc_not_cpg_island: list of lists with the start point and end point of all the not cpg island
    :return: two lists: a list of the seqs and a list of their labels (accordingly)
    """
    seqs_list = list()
    labels_list = list()
    for not_cpg in loc_not_cpg_island:
        seqs_list.append(chrom[not_cpg[0]:not_cpg[1]].upper())
        labels_list.append(0)

    for cpg in loc_cpg_island:
        if cpg[1] - cpg[0] < THRESHOLD:
            continue
        seqs_list.append(chrom[cpg[0]:cpg[0] + SEQ_LEN].upper())
        labels_list.append(1)

    return seqs_list, labels_list


if __name__ == '__main__':
    """ sys.argv[1] = the wanted chromosome """
    loc_cpg_island, loc_not_cpg_island = get_locations()
    chrom = get_chrom(sys.argv[1])
    seqs_list, labels_list = set_data(chrom, loc_cpg_island, loc_not_cpg_island)

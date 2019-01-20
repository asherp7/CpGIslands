import os
import pickle
from utils import *
import sys

NUM_OF_CHROM = 6
THRESHOLD = 500
SEQ_LEN = 1001


def get_locations(genome_tag_range, chrom_name):
    """
    :return: Two lists of lists with the location of cpg island and not
    """
    # d1 = genome_tag_range.get_untagged_ranges('chr1', distance=1000, range_len=2000)
    # d1 = genome_tag_range.get_random_untagged_ranges('chr1', ranges_num=100, range_len=10000)
    d1 = genome_tag_range.get_tagged_ranges(chrom_name)
    loc_cpg_island = list()
    for start, end in zip(d1[RANGE_START], d1[RANGE_END]):
        # print('%d: %d -- %d' % (count, start, end))
        loc_cpg_island.append([start, end])
    loc_not_cpg_island = list()
    d1 = genome_tag_range.get_untagged_ranges(chrom_name, THRESHOLD, SEQ_LEN)
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
    # :param chrom_name: the name of the wanted chorom
    # :param chrom_string: the wanted chromosome as string
    # :param loc_cpg_island: locations of the cpg islands
    # :param loc_not_cpg_island: locations of the not cpg islands
    # :param genome_tag_range: object of the class utils
    # :return: for each chromosom, a tuple of three which include : the Str of the seq, his label, and its nucleotide
        label for each letter
    # """

    tuple_of_data = list()
    for not_cpg in loc_not_cpg_island:
        a = chrom_string[not_cpg[0]:not_cpg[1]].upper()
        b = "0"
        start, length = not_cpg[0], not_cpg[1] - not_cpg[0]
        c = genome_tag_range.get_tag_for_range(chrom_name, start, length)
        tuple_of_data.append((a, b, c))
    for cpg in loc_cpg_island:
        if cpg[1] - cpg[0] < THRESHOLD:
            continue
        a = chrom_string[cpg[0]:cpg[0] + SEQ_LEN].upper()
        b = "1"
        start, length = cpg[0], cpg[1] - cpg[0]
        c = genome_tag_range.get_tag_for_range(chrom_name, start, length)
        tuple_of_data.append((a, b, c))
    return tuple_of_data


if __name__ == '__main__':
    """ sys.argv[1] = directory name for the split data folder
        sys.argv[2] = the location of the CGI.hg19.bed file
        sys.argv[3] = the location of the hg19.chrom.sizes file"""

    dir_name = sys.argv[1]
    genome_tag_range = GenomeTagRange(sys.argv[2], sys.argv[3])
    all_data = list()

    i = 0
    for chrom_name in os.listdir(dir_name):
        if chrom_name != "chr1" and i < NUM_OF_CHROM:
            print(chrom_name)
            chrom_string = get_chrom(chrom_name)
            loc_cpg_island, loc_not_cpg_island = get_locations(genome_tag_range, chrom_name)
            one_chrom = set_data(chrom_name, chrom_string, loc_cpg_island, loc_not_cpg_island, genome_tag_range)
            all_data.extend(one_chrom)
            i += 1

    with open('all_data', 'wb') as f:
        pickle.dump(all_data, f)


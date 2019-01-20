from utils import *
from Bio import SeqIO, SeqUtils
import numpy as np


if __name__ == "__main__":
    file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\test.txt'

    new_file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\test_new.txt'
    file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.fa'
    new_file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.fa'

    genome_tag_range = GenomeTagRange(r'C:\Users\User\Documents\Thirdyear\Tommy_course\CpGIslands\hack_proj\data\CGI.hg19.bed',
                                      r'C:\Users\User\Documents\Thirdyear\Tommy_course\CpGIslands\hack_proj\data\hg19.chrom.sizes')
    # d1 = genome_tag_range.get_untagged_ranges('chr1', distance=1000, range_len=2000)
    d1 = genome_tag_range.get_random_untagged_ranges('chr1', ranges_num=100, range_len=10000)
    # d1 = genome_tag_range.get_tagged_ranges('chr1')
    d1 = genome_tag_range.get_untagged_ranges('chr1', 500, 1001)
    count = 1
    for start, end in zip(d1[RANGE_START], d1[RANGE_END]):
        print('%d: %d -- %d' % (count, start, end))
        count += 1

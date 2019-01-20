from utils import *
from Bio import SeqIO, SeqUtils
import numpy as np


if __name__ == "__main__":
    file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\test.txt'

    new_file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\test_new.txt'
    file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.fa'
    new_file_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.fa'

    a = GenomeTagRange(r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\CGI.hg19.bed')
    print(a.get_names())
    ret = a.get_tag_for_range('chrUn_gl000225', 30300, 3)
    # ret = a.get_tag_for_range('chr1', 805199, 1451)

    # ret = a.get_tag_for_range('chr1', 29809, 3)
    ret = np.array(ret)
    print(ret)
    print(np.sum(ret))
    # with open(file_path, 'rt') as fd:
    #     for record in SeqIO.parse(fd, 'fasta'):
    #         with open(new_file_path + '.' + record.id, 'wt') as new_fd:
    #             new_fd.write(record.seq._data)
    #         print(record.id)
    #
    #     SeqIO.read(fd, 'fasta')
    #     start = 0
    #     found_chr_2 = False
    #     counter = 1
    #     while not found_chr_2:
    #         print('%d: Read %d to %d' % (counter, start, start + 5000))
    #         s = read_range(fd, start, 5000)
    #         if 'chr3' in s:
    #             found_chr_2 = True
    #         start += 5000
    #         counter += 1
    #     print(s)

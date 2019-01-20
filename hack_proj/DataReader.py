import pandas as pd
from Bio import SeqIO
import os


class DataReader:
    def __init__(self, data_path):
        self.data_path = data_path
        self.split_data = data_path + 'split_data/'
        self.chrome_sizes = pd.read_csv(data_path + 'hg19.chrom.sizes', sep='\t', names=['CHROMOSOME', 'SIZE'])
        self.cgi_bed = pd.read_csv(data_path + 'CGI.hg19.bed', sep='\t', names=['CHROMOSOME', 'START', 'END', 'SCORE'])
        self.hg19fa_path = data_path + 'hg19.fa'
        self.hg19_cpg_bed_path = data_path + 'hg19.CpG.bed'
        self.chromosome_list = []
        self.data_file_object = None

        # init functions:
        self._add_cpg_island_lengths()
        self._set_chromosome_list()
        self.open_data_file_object()

    def __del__(self):
        if self.data_file_object:
            self.data_file_object.close()

    def open_data_file_object(self):
        self.data_file_object = open(self.hg19fa_path)

    def _set_chromosome_list(self):
        self.chromosome_list = self.chrome_sizes['CHROMOSOME'].tolist()

    def _add_cpg_island_lengths(self):
        self.cgi_bed['LENGTH'] = self.cgi_bed['END'] - self.cgi_bed['START']

    def get_sequence(self, chromosome, start_index, n):
        with open(self.split_data + chromosome, 'rt') as data_file:
            data_file.seek(start_index)
            return data_file.read(n)

    # def read_fasta(self):
    #     record_list = list(SeqIO.parse(self.hg19fa_path, "fasta"))
    #     print(record_list[0], len(record_list[0]))

    # def readInChunks(self, chunkSize=1024):
    #     """
    #     Lazy function to read a file piece by piece.
    #     Default chunk size: 1kB.
    #     """
    #     while True:
    #         data = self.data_file_object.read(chunkSize)
    #         if not data:
    #             break
    #         yield data

    # def devideToChunks(self, chunkSize=1024):
    #     """
    #     function that splits large file into small files of size chuncksize
    #     :param chunkSize: the number of chars each file will contain
    #     """
    #
    #     for idx, chunk in enumerate(self.readInChunks(chunkSize)):
    #         with open("Output.txt", "w") as text_file:
    #             text_file.write(chunk)

    def devideDataToChromosomes(self):
        """
        split large data_file into one file per chromosome
        :param new_file_path: path to write the splitted data
        :return:
        """
        if not os.path.exists(self.split_data):
            os.mkdir(self.split_data)
        with open(self.hg19fa_path, 'rt') as fd:
            for record in SeqIO.parse(fd, 'fasta'):
                with open(self.split_data + record.id, 'wt') as new_fd:
                    new_fd.write(record.seq._data)
                print(record.id)



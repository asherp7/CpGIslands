from DataReader import DataReader

data_reader = DataReader('data/')
# To divide data into seperate files for each chromosome:
data_reader.devideDataToChromosomes()

# data_reader.read_fasta()
# my_gen = data_reader.readInChunks(128)
# print(next(my_gen))
print('Number of chromosomes: ', len(data_reader.chromosome_list))
print('AVG chromosome size: ', int(data_reader.chrome_sizes['SIZE'].mean()))
print('Number of cpg islands: ', len(data_reader.cgi_bed.index))
print('AVG CpG island size: ', int(data_reader.cgi_bed['LENGTH'].mean()))
print('hg19 size (93 x AVG_CHROME_SIZE): ', len(data_reader.chromosome_list) * int(data_reader.chrome_sizes['SIZE'].mean()))
print('chromosome list:', data_reader.chromosome_list)
sequence = data_reader.get_sequence('chr1', 10000, 100)
print('get sequence example: ', sequence)
print('histogram of sequence:', data_reader.get_k_mer_histogram(sequence, 5, True))




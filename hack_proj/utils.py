import csv
import bisect
import numpy as np


RANGE_START = 'start'
RANGE_END = 'end'
RANGE_LEN = 'length'
RANGE_SCORE = 'score'


def read_range(fd, start, length):
    fd.seek(start)
    ret = fd.read(length)
    return ret


class GenomeTagRange(object):
    def __init__(self, tag_file_path, chromosome_size_file_path):
        self.tag_file_path = tag_file_path
        self.chromosome_size_file_path = chromosome_size_file_path
        self.data_dict = None
        self.sizes_dict = None
        self._init_data()

    def _init_data(self):
        self.data_dict = {}
        with open(self.tag_file_path) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')

            for row in tsv_reader:
                if row[0] not in self.data_dict:
                    self.data_dict[row[0]] = {RANGE_START: [], RANGE_END: [], RANGE_SCORE: []}
                self.data_dict[row[0]][RANGE_START].append(int(row[1]))
                self.data_dict[row[0]][RANGE_END].append(int(row[2]))
                self.data_dict[row[0]][RANGE_SCORE].append(int(row[3]))
        self.sizes_dict = self._init_sizes(self.chromosome_size_file_path)

    def _init_sizes(self, file_path):
        sizes_dict = {}
        with open(file_path) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')
            for row in tsv_reader:
                sizes_dict[row[0]] = int(row[1])
        return sizes_dict

    def get_names(self):
        return list(self.data_dict.keys())

    def get_tagged_ranges(self, chr):
        return self.data_dict[chr]

    def get_distance_from_tagged(self, chr, start, end):
        i = bisect.bisect_right(self.data_dict[chr]['start'], start) - 1
        if i < 0:
            if end < self.data_dict[chr]['start'][0]:
                return self.data_dict[chr]['start'][0] - end
            else:
                return 0
        else:
            if start < self.data_dict[chr]['end'][i]:
                return 0

            if i < (len(self.data_dict[chr]['start']) - 1) and end >= self.data_dict[chr]['start'][i+1]:
                return 0

            if i == (len(self.data_dict[chr]['start']) - 1):
                return start - self.data_dict[chr]['end'][i]
            else:
                return min(start - self.data_dict[chr]['end'][i], self.data_dict[chr]['start'][i+1] - end)

    def get_untagged_ranges(self, chr, distance=1000, range_len=None):
        new_dict = {RANGE_START: [], RANGE_END: [], RANGE_LEN: []}

        for i in range(len(self.data_dict[chr][RANGE_START])):
            if range_len is None:
                curr_len = self.data_dict[chr][RANGE_END][i] - self.data_dict[chr][RANGE_START][i]
            else:
                curr_len = range_len
            curr_start = self.data_dict[chr][RANGE_START][i] - curr_len - distance
            curr_end = curr_start + curr_len
            if self.get_distance_from_tagged(chr, curr_start, curr_end) >= distance:
                new_dict[RANGE_START].append(curr_start)
                new_dict[RANGE_END].append(curr_start + curr_len)
                new_dict[RANGE_LEN].append(curr_len)

            curr_start = self.data_dict[chr][RANGE_END][i] + distance
            curr_end = curr_start + curr_len
            if self.get_distance_from_tagged(chr, curr_start, curr_end) >= distance:
                new_dict[RANGE_START].append(curr_start)
                new_dict[RANGE_END].append(curr_start + curr_len)
                new_dict[RANGE_LEN].append(curr_len)

        return new_dict

    def get_random_untagged_ranges(self, chr, ranges_num, range_len, distance=1000):
        new_dict = {RANGE_START: [], RANGE_END: [], RANGE_LEN: []}
        tries = 0
        max_tries = ranges_num * ranges_num
        curr_size = self.sizes_dict[chr]
        min_base = int(0.00001 * curr_size)
        max_base = curr_size - range_len - 1
        while len(new_dict[RANGE_START]) < ranges_num and tries < max_tries:
            tries += 1
            curr_start = np.random.randint(min_base, max_base)
            curr_end = curr_start + range_len
            if self.get_distance_from_tagged(chr, curr_start, curr_end) >= distance:
                new_dict[RANGE_START].append(curr_start)
                new_dict[RANGE_END].append(curr_end)
                new_dict[RANGE_LEN].append(range_len)
        return new_dict

    def get_tag_for_range(self, chr, start, length):
        if length <= 0:
            return []

        i = bisect.bisect_right(self.data_dict[chr]['start'], start) - 1
        # Not in cpg islands range
        if i < 0:
            curr_island_tag = [0] * (self.data_dict[chr]['start'][0] - start)
            extension_list = self.get_tag_for_range(chr, self.data_dict[chr]['start'][0], length - (self.data_dict[chr]['start'][0] - start))
            curr_island_tag.extend(extension_list)
            return curr_island_tag

        range_end = start + length
        curr_start = self.data_dict[chr]['start'][i]
        curr_end = self.data_dict[chr]['end'][i]

        # Starts after the current tagged range
        if start > curr_end:
            if i == len(self.data_dict[chr]['start']) - 1:
                return [0] * length
            else:
                if range_end < self.data_dict[chr]['start'][i+1]:
                    return [0] * length
                else:
                    curr_island_tag = [0] * (self.data_dict[chr]['start'][i+1] - start)
                    extension_list = self.get_tag_for_range(chr, self.data_dict[chr]['start'][i+1], length - (self.data_dict[chr]['start'][i+1] - start))
                    curr_island_tag.extend(extension_list)
                    return curr_island_tag

        # Start at current end
        if start == curr_end:
            curr_island_tag = [0]
            extension_list = self.get_tag_for_range(chr, start + 1, length - 1)
            curr_island_tag.extend(extension_list)
            return curr_island_tag

        # Finish before current range ends
        if curr_end >= range_end:
            return [1] * length

        # Finish after current end
        if curr_end < range_end:
            curr_island_tag = [1] * (curr_end - start)
            extension_list = self.get_tag_for_range(chr, curr_end, length - (curr_end - start))
            curr_island_tag.extend(extension_list)
            return curr_island_tag

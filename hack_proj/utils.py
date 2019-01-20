import csv
import bisect


def read_range(fd, start, length):
    fd.seek(start)
    ret = fd.read(length)
    return ret


class GenomeTagRange(object):
    def __init__(self, tag_file_path):
        self.tag_file_path = tag_file_path
        self.data_dict = None
        self._init_data()

    def _init_data(self):
        self.data_dict = {}
        with open(self.tag_file_path) as tsv_file:
            tsv_reader = csv.reader(tsv_file, delimiter='\t')

            for row in tsv_reader:
                if row[0] not in self.data_dict:
                    self.data_dict[row[0]] = {'start': [], 'end': [], 'score': []}
                self.data_dict[row[0]]['start'].append(int(row[1]))
                self.data_dict[row[0]]['end'].append(int(row[2]))
                self.data_dict[row[0]]['score'].append(int(row[3]))

    def get_names(self):
        return list(self.data_dict.keys())

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

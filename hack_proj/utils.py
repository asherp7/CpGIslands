import csv
import bisect
import numpy as np
from sklearn.metrics import roc_curve, auc, classification_report
import matplotlib.pyplot as plt
from sklearn.svm import SVC
import collections

import DataReader

RANGE_START = 'start'
RANGE_END = 'end'
RANGE_LEN = 'length'
RANGE_SCORE = 'score'

ALL_BASES = ['A', 'C', 'G', 'T']
UNKNOWN_BASE = 'N'


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


class DataExtractor(object):
    def __init__(self, data_reader, genome_tag_reader):
        self.data_reader = data_reader  # type: DataReader.DataReader
        self.genome_tag_reader = genome_tag_reader  # type: GenomeTagRange

    def get_data_from_ranges(self, chr, ranges_dict, pad_size=0, with_unknown=False):
        all_seq = []
        all_labels = []
        seq_num = len(ranges_dict[RANGE_START])
        for i in range(seq_num):
            curr_start = ranges_dict[RANGE_START][i] - pad_size
            curr_end = ranges_dict[RANGE_END][i] + pad_size

            curr_labels = self.genome_tag_reader.get_tag_for_range(chr, curr_start, curr_end - curr_start)
            curr_seq = self.data_reader.get_sequence('hg19.fa.' + chr, curr_start, curr_end - curr_start)

            if with_unknown or curr_seq.find(UNKNOWN_BASE) < 0:
                all_seq.append(curr_seq)
                all_labels.append(curr_labels)
            else:
                print('%d-%s: Found %s in sequence at %d' % ((i+1), chr, UNKNOWN_BASE, curr_seq.find(UNKNOWN_BASE) + curr_start))
        if len(all_seq) < seq_num:
            print('Returned %d sequqnces instead of %d. Missed %d sequences.' % (len(all_seq), seq_num, seq_num - len(all_seq)))
        return all_seq, all_labels

    def get_base_tagged_seq(self, chr, seq_num=1, pad_size=1000, with_unknown=False):
        all_ranges = self.genome_tag_reader.get_tagged_ranges(chr)
        seq_num = min(seq_num, len(all_ranges[RANGE_START]))
        all_ranges[RANGE_START] = all_ranges[RANGE_START][:seq_num]
        all_ranges[RANGE_END] = all_ranges[RANGE_END][:seq_num]
        return self.get_data_from_ranges(chr, all_ranges, pad_size, with_unknown)

    def get_non_island_random_data(self, chr, seq_num=1, seq_len=1000, min_distance=1000, with_unknown=False):
        all_ranges = self.genome_tag_reader.get_random_untagged_ranges(chr, seq_num, seq_len, distance=min_distance)
        return self.get_data_from_ranges(chr, all_ranges, 0, with_unknown)

    def get_near_island_data(self, chr, seq_num=1, seq_len=1000, distance=1000, with_unknown=False):
        all_ranges = self.genome_tag_reader.get_untagged_ranges(chr, distance, seq_len)
        seq_num = min(seq_num, 2 * len(all_ranges[RANGE_START]))
        all_ranges[RANGE_START] = all_ranges[RANGE_START][:seq_num]
        all_ranges[RANGE_END] = all_ranges[RANGE_END][:seq_num]
        return self.get_data_from_ranges(chr, all_ranges, 0, with_unknown)

def get_results_from_dict_list(dict_list, label_key='label', pred_key='score'):
    y = []
    y_pred = []
    for d in dict_list:
        y.append(d[label_key])
        y_pred.append(d[pred_key])
    return y, y_pred


def show_roc_curve(y, y_score, title='ROC Curve', show_grid=True, print_data=True):
    fpr, tpr, thresholds = roc_curve(y, y_score)
    roc_auc = auc(fpr, tpr)

    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc="lower right")
    plt.grid(show_grid)

    # create the axis of thresholds (scores)
    ax2 = plt.gca().twinx()
    ax2.plot(fpr, thresholds, markeredgecolor='r', linestyle='dashed', color='g')
    ax2.set_ylabel('Threshold', color='g')
    ax2.set_ylim([thresholds[-1], thresholds[0]])
    ax2.set_xlim([fpr[0], fpr[-1]])
    if print_data:
        for i in range(fpr.size):
            print('FPR: %.3f, TPR: %.3f: Threshold = %.4f' % (fpr[i], tpr[i], thresholds[i]))
    plt.show()


def get_prediction_stats(y, y_pred, target_names=('Regular', 'CPG')):
    return classification_report(y, y_pred, target_names=target_names)


def print_prediction_stats(y, y_pred, target_names=('Regular', 'CPG')):
    print(get_prediction_stats(y, y_pred, target_names))


def get_seq_windows(seq, window_size, base_labels=None):
    all_seq = []
    all_labels = []
    side = int(window_size / 2)
    for i in range(side, len(seq) - side):
        all_seq.append(seq[i-side:i+side])
        if base_labels is not None:
            all_labels.append(base_labels[i-side:i+side])
    new_seq = seq[side:len(seq) - side]
    if base_labels is None:
        return new_seq, all_seq
    else:
        new_lables = base_labels[side:len(seq) - side]
        return new_seq, all_seq, new_lables, all_labels


def get_random_sample(all_seq, all_labels=None, ratio=0.01, sample_size=None):
    if sample_size is None:
        sample_size = max(int(len(all_seq) * ratio), 1)
    idx = np.random.choice(len(all_seq), sample_size, replace=False)
    new_seq = [all_seq[j] for j in idx]
    if all_labels is None:
        return all_seq
    else:
        new_labels = [all_labels[j] for j in idx]
        return all_seq, all_labels


def count_substr(s, sub_str):
    count = 0
    start = 0
    while True:
        start = s.find(sub_str, start) + 1
        if start > 0:
            count += 1
        else:
            return count


class IslandClassifier(object):
    def __init__(self, threshold):
        self.threshold = threshold

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_seq_prob(self, all_seq):
        pass

    def get_seq_classification(self, all_seq):
        pass


class SubstrFreqClassifier(IslandClassifier):
    def __init__(self, substr_list, seq_window_len, threshold=None):
        IslandClassifier.__init__(self, threshold)
        self.seq_window_len = seq_window_len
        self.svc_clf = None  # type: SVC
        self.substr_list = substr_list

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_seq_freqs(self, seq):
        freqs = []

        for i, s in enumerate(self.substr_list):
            if len(s) == 1:
                freqs.append(seq.count(s) / len(seq))
            else:
                freqs.append(count_substr(seq, s) / (len(seq) - len(s) + 1))

        return freqs

    def lst2arr(self, lst):
        arr = np.array(lst)
        if len(arr.shape) == 1:
            arr = np.reshape(arr, (-1, 1))
        return arr

    def get_data_from_seq_list(self, all_seq):
        x = []
        for i in range(len(all_seq)):
            seq = all_seq[i]
            freqs = self.get_seq_freqs(seq)
            x.append(freqs)

        return x

    def fit(self, all_seq, y, use_window=True, seq_ratio=0.01):
        x = []
        new_y = []
        print('Create data to fit...')
        for i in range(len(all_seq)):
            seq = all_seq[i]
            if not use_window:
                freqs = self.get_seq_freqs(seq)
                x.append(freqs)
            else:
                new_seq, all_windows_seq = get_seq_windows(seq, self.seq_window_len)
                all_windows_seq = get_random_sample(all_windows_seq, all_labels=None, ratio=seq_ratio)
                curr_x = self.get_data_from_seq_list(all_windows_seq)
                x.extend(curr_x)
                new_y.extend([y[i]] * len(curr_x))

        if use_window:
            y = new_y

        x = self.lst2arr(x)
        self.svc_clf = SVC(kernel='linear', probability=True, max_iter=10000)
        print('Fit data (%d samples)...' % (x.shape[0]))
        self.svc_clf.fit(x, np.array(y))

    def get_seq_prob(self, all_seq):
        x = self.get_data_from_seq_list(all_seq)
        y_pred = self.svc_clf.predict_proba(self.lst2arr(x))
        return y_pred[:, 1]

    def classify_seq(self, seq, labels=None):
        new_seq, all_seq, new_lables, all_labels = get_seq_windows(seq, window_size=self.seq_window_len,
                                                                   base_labels=labels)
        y_prob = self.get_seq_prob(all_seq)
        y_pred = np.zeros_like(y_prob)
        y_pred[y_prob >= self.threshold] = 1
        return y_pred


class CGFreqClassifier(IslandClassifier):
    def __init__(self, seq_window_len, cg_threshold=None):
        IslandClassifier.__init__(self, cg_threshold)
        self.seq_window_len = seq_window_len
        self.svc_clf = None  # type: SVC

    def set_threshold(self, threshold):
        self.threshold = threshold

    def get_seq_counts(self, seq):
        count_dict = {}
        for c in ALL_BASES:
            count_dict[c] = 0

        for i in range(len(seq)):
            try:
                count_dict[seq[i]] += 1
            except:
                continue
        return count_dict

    def get_data_from_seq_list(self, all_seq):
        x = []
        for i in range(len(all_seq)):
            seq = all_seq[i]
            count_dict = self.get_seq_counts(seq)
            x.append(((count_dict.get('G', 0) + count_dict.get('C', 0)) / len(seq)))

        x = np.reshape(np.array(x), (-1, 1))
        return x

    def fit(self, all_seq, y):
        x = []
        new_y = []
        print('Create data to fit...')
        for i in range(len(all_seq)):
            seq = all_seq[i]
            new_seq, all_windows_seq = get_seq_windows(seq, self.seq_window_len)
            idx = np.random.choice(len(all_windows_seq), max(int(len(all_windows_seq) * 0.01), 1), replace=False)
            all_windows_seq = [all_windows_seq[j] for j in idx]
            curr_x = self.get_data_from_seq_list(all_windows_seq)
            # curr_x = np.random.choice(np.squeeze(curr_x), max(int(curr_x.size * 0.01), 1), replace=False)
            for val in curr_x:
                x.append(val)
                new_y.append(y[i])
            # for win_seq in all_windows_seq:
            #     count_dict = self.get_seq_counts(win_seq)
            #     x.append(((count_dict.get('G', 0) + count_dict.get('C', 0)) / len(win_seq)))
            #     new_y.append(y[i])

        x = np.reshape(np.array(x), (-1, 1))
        self.svc_clf = SVC(kernel='linear', probability=True, max_iter=10000)
        print('Fit data...')
        self.svc_clf.fit(np.array(x), np.array(new_y))

    def get_seq_prob(self, all_seq):
        x = self.get_data_from_seq_list(all_seq)
        y_pred = self.svc_clf.predict_proba(x)
        return y_pred[:, 1]

    def classify_seq(self, seq, labels=None):
        new_seq, all_seq, new_lables, all_labels = get_seq_windows(seq, window_size=self.seq_window_len,
                                                                   base_labels=labels)
        y_prob = self.get_seq_prob(all_seq)
        y_pred = np.zeros_like(y_prob)
        y_pred[y_prob >= self.threshold] = 1
        return y_pred

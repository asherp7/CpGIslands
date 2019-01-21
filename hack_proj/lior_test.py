from utils import *
from Bio import SeqIO, SeqUtils
import numpy as np
import pickle
import time
from markov_model import *


def create_data_01(train_list, test_list):
    data_dir_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data/'
    sizes_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.chrom.sizes'
    tags_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\CGI.hg19.bed'
    output_dir = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01'

    genome_tag_range = GenomeTagRange(tags_path, sizes_path)
    reader = DataReader.DataReader(data_dir_path)
    data_extractor = DataExtractor(reader, genome_tag_range)

    print('Create train_islands...')
    pkl_path = os.path.join(output_dir, 'train_islands.pkl')
    if not os.path.isfile(pkl_path):
        train_islands, _ = create_island_overlap_data(train_list, data_extractor, pad_size=0, seq_num_list=None)
        pickle.dump(train_islands, open(os.path.join(output_dir, 'train_islands.pkl'), 'wb'))
        train_islands.clear()

    print('Create train_random...')
    pkl_path = os.path.join(output_dir, 'train_random.pkl')
    if not os.path.isfile(pkl_path):
        train_random, _ = create_non_island_random_data(train_list, data_extractor, seq_len=2001, min_distance=5000, seq_num_list=500)
        pickle.dump(train_random, open(os.path.join(output_dir, 'train_random.pkl'), 'wb'))
        train_random.clear()

    print('Create train_near...')
    pkl_path = os.path.join(output_dir, 'train_near.pkl')
    if not os.path.isfile(pkl_path):
        train_near, _ = create_near_island_data(train_list, data_extractor, seq_num_list=500, distance=1000, seq_len=2001)
        pickle.dump(train_near, open(os.path.join(output_dir, 'train_near.pkl'), 'wb'))
        train_near.clear()

    print('Create test padded_islands')
    pkl_path = os.path.join(output_dir, 'test_padded_islands_02.pkl')
    if not os.path.isfile(pkl_path):
        test_padded_islands, labels = create_island_overlap_data(test_list, data_extractor, pad_size=1500, seq_num_list=None)
        full_str = ''
        full_lables = []
        for i in range(len(test_padded_islands)):
            full_str = full_str + test_padded_islands[i]
            full_lables.extend(labels[i])
        pickle.dump((full_str, full_lables), open(pkl_path, 'wb'))


def test_new_markov():
    markov_order = 2
    win_size = 51
    island_markov_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\markov_model_island_%d.pkl' % markov_order
    other_markov_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\markov_model_other_%d.pkl' % markov_order
    test_seq_pkl_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\chr1_125124_145563.pkl'
    test_seq_pkl_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\test_padded_islands_02.pkl'

    island_markov = pickle.load(open(island_markov_path, 'rb')) # type: MarkovModel
    llr_markov = pickle.load(open(island_markov_path, 'rb')) # type: MarkovModel
    other_markov = pickle.load(open(other_markov_path, 'rb')) # type: MarkovModel
    seq, labels = pickle.load(open(test_seq_pkl_path, 'rb'))
    seq = seq[:100000]
    labels = labels[:100000]
    llr_markov.log_prob_mat = island_markov.log_prob_mat - other_markov.log_prob_mat
    llr_markov.log_transition_mat = island_markov.log_transition_mat - other_markov.log_transition_mat
    print_transition_mat(island_markov)

    llr_score = llr_markov.get_ll([seq])
    llr_score_island = island_markov.get_ll([seq])
    llr_score_other = other_markov.get_ll([seq])

    print('N Count: %d' % count_substr(seq, 'N'))
    new_labels = labels[llr_markov.order:]
    window_score, window_labels = apply_window(llr_score, new_labels, win_size=win_size)

    show_roc_curve(window_labels, window_score, print_data=False)
    y_pred = np.zeros_like(window_score)
    y_pred[window_score > 17.3] = 1
    print_prediction_stats(window_labels, y_pred)
    idx = np.array(range(window_score.size))

    plt.figure()
    plot_scores_and_labels(window_score, window_labels, 'LLR', 'Log Likeklihood Ratio - Markov order = %d, Win size = %d' % (markov_order, win_size))
    # plt.figure()
    # window_score_other, window_labels_other = apply_window(llr_score_other, new_labels, win_size=win_size)
    # plot_scores_and_labels(window_score_other, window_labels, 'Other LL', 'Other Sequence Log Likelihood - Markov order %d, Win size = %d' % (markov_order, win_size))
    # plt.figure()
    # window_score_island, window_labels_island = apply_window(llr_score_island, new_labels, win_size=win_size)
    # plot_scores_and_labels(window_score_island, window_labels, 'Island LL', 'Island Sequence Log Likelihood - Markov order %d, Win size = %d' % (markov_order, win_size))
    plt.show()


def plot_scores_and_labels(score, labels, score_type, title):
    idx = np.array(range(score.size))
    # plt.scatter(range(score.size), normalize_arr(score), marker='.', label=score_type)
    plt.plot(range(score.size), normalize_arr(score), label=score_type)
    plt.scatter(idx[labels == 1], labels[labels == 1], marker='.', label='CPG Island', c='orange')
    plt.scatter(idx[labels == 0], labels[labels == 0], marker='.', label='Other', c='green')
    plt.xlabel('DNA Sequence location')
    plt.ylabel('CPG Island score')
    plt.title(title)
    plt.legend()


def train_markov_models():
    island_paths = [r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\train_islands.pkl']
    other_paths = [r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\train_random.pkl',
                   r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01\train_near.pkl']
    output_dir = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\create_data_01'

    markov_order = 5
    print('----------------Island data--------------')
    island_data = []
    for p in island_paths:
        island_data.extend(pickle.load(open(p, 'rb')))
    markov_model_island = MarkovModel(ALL_BASES, order=markov_order)
    # markov_model_island.fit_transition(island_data[:10])
    # print_transition_mat(markov_model_island)
    # markov_model_island.fit_transition(island_data[:100])
    # print_transition_mat(markov_model_island)
    # markov_model_island.fit_transition(island_data[:1000])
    # print_transition_mat(markov_model_island)
    markov_model_island.fit_transition(island_data)
    print_transition_mat(markov_model_island)
    pickle.dump(markov_model_island, open(os.path.join(output_dir, 'markov_model_island_%d.pkl' % markov_order), 'wb'))
    print('----------------Other data--------------')
    other_data = []
    for p in other_paths:
        other_data.extend(pickle.load(open(p, 'rb')))
    markov_model_other = MarkovModel(ALL_BASES, order=markov_order)
    # markov_model_other.fit_transition(other_data[:10])
    # print_transition_mat(markov_model_other)
    # markov_model_other.fit_transition(other_data[:100])
    # print_transition_mat(markov_model_other)
    # markov_model_other.fit_transition(other_data[:1000])
    # print_transition_mat(markov_model_other)
    markov_model_other.fit_transition(other_data)
    print_transition_mat(markov_model_other)
    pickle.dump(markov_model_other, open(os.path.join(output_dir, 'markov_model_other_%d.pkl' % markov_order), 'wb'))


def extract_data_01():
    data_dir_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data/'
    sizes_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.chrom.sizes'
    tags_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\CGI.hg19.bed'
    data_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\clf data\all_data.pkl'

    genome_tag_range = GenomeTagRange(tags_path, sizes_path)
    reader = DataReader.DataReader(data_dir_path)
    data_extractor = DataExtractor(reader, genome_tag_range)

    start = 125124
    end = 145563
    chr = 'chr1'
    ranges_dict = {RANGE_START: [start], RANGE_END: [end]}

    pkl_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\tmp\%s_%d_%d.pkl' % (chr, start, end)
    all_seq, all_labels = data_extractor.get_data_from_ranges(chr, ranges_dict, pad_size=0, with_unknown=True)
    pickle.dump((all_seq[0], all_labels[0]), open(pkl_path, 'wb'))


def test_all_01():
    data_dir_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data/'
    sizes_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.chrom.sizes'
    tags_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\CGI.hg19.bed'
    data_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\clf data\all_data.pkl'

    genome_tag_range = GenomeTagRange(tags_path, sizes_path)
    reader = DataReader.DataReader(data_dir_path)
    data_extractor = DataExtractor(reader, genome_tag_range)

    t = time.time()
    all_seq, all_labels = create_island_overlap_data(train_chrs, data_extractor, pad_size=1000, seq_num_list=None)
    print('Found %d sequences' % len(all_seq))
    print('create_island_overlap_data Time: %f seconds' % (time.time() - t))

    t = time.time()
    all_seq, all_labels = create_near_island_data(train_chrs, data_extractor, None)
    print('Found %d sequences' % len(all_seq))
    print('create_near_island_data Time: %f seconds' % (time.time() - t))


    t = time.time()
    all_seq, all_labels = create_non_island_random_data(train_chrs, data_extractor, 2000)
    print('Found %d sequences' % len(all_seq))
    print('create_non_island_random_data Time: %f seconds' % (time.time() - t))

    exit(1)
    all_data = pickle.load(open(data_path, 'rb'))
    all_seq = []
    all_label = []
    for seq, label, seq_label in all_data:
        all_seq.append(seq)
        all_label.append(int(label))
    # clf = CGFreqClassifier(seq_window_len=501)
    clf = SubstrFreqClassifier(substr_list=['G', 'C'], seq_window_len=501)
    clf = SubstrFreqClassifier(substr_list=['GC', 'CG'], seq_window_len=501)
    train_ratio = 0.9
    test_ratio = 0.9
    clf.fit(all_seq[:int(len(all_seq) * train_ratio)], all_label[:int(len(all_seq) * train_ratio)])
    print('Predict test...')

    y = all_label[int(len(all_seq) * test_ratio):]
    test_data = all_seq[int(len(all_seq) * test_ratio):]
    y_pred = clf.get_seq_prob(test_data)
    show_roc_curve(y, y_pred)

    thresh = 0.3
    clf.set_threshold(thresh)

    # all_seq, all_labels = extractor.get_base_tagged_seq('chr1', pad_size=2)


def test_data_01():
    data_dir_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data/'
    sizes_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\hg19.chrom.sizes'
    tags_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\CGI.hg19.bed'

    genome_tag_range = GenomeTagRange(tags_path, sizes_path)
    reader = DataReader.DataReader(data_dir_path)
    extractor = DataExtractor(reader, genome_tag_range)
    all_seq, all_labels = extractor.get_base_tagged_seq('chr1', pad_size=2)

    print(all_seq)
    print(all_labels)


def test_classifier_01():
    data_path = r'C:\liorz\school\76558 Algorithms in Computational Biology\hackathon\data\clf data\all_data.pkl'
    all_data = pickle.load(open(data_path, 'rb'))
    all_seq = []
    all_label = []
    for seq, label, seq_label in all_data:
        all_seq.append(seq)
        all_label.append(int(label))
    clf = CGFreqClassifier(100)
    clf.fit(all_seq[:int(len(all_seq) * 0.9)], all_label[:int(len(all_seq) * 0.9)])
    y = all_label[int(len(all_seq) * 0.9):]
    test_data = all_seq[int(len(all_seq) * 0.9):]
    y_pred = clf.get_seq_prob(test_data)
    show_roc_curve(y, y_pred)


if __name__ == "__main__":
    # train_markov_models()
    # create_data_01(train_chrs, test_chrs)
    test_new_markov()
    # extract_data_01()
    exit(1)

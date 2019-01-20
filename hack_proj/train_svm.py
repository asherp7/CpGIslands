import itertools
from collections import Counter
import numpy as np
np.set_printoptions(threshold=np.nan)
from sklearn import svm
import pickle

CHROMOSOME = 'chr1'
K_MER_LEN = 6
TESTING_RATIO = 0.1


def svm_testing_loss(clf, test_data, test_labels):
    predictions = clf.predict(test_data)
    fn, fp, tn, tp = 0, 0, 0, 0
    for predction, label in zip(predictions, test_labels):
        if predction == label == 0:
            tn += 1
        elif predction == label == 1:
            tp += 1
        elif predction == 1 and label == 0:
            fp += 1
        elif predction == 0 and label == 1:
            fn += 1
    recall = tp / (tp + fn)
    precision = tp / (tp + fp)
    num_predections = len(predictions)

    # print svm results:
    print('***********************')
    print('*                     *')
    print('*     SVM RESULTS:    *')
    print('*                     *')
    print('* RECALL = ', recall)
    print('* PRECISION = ', precision)
    print('* TP = ', tp)
    print('* TN = ', tn)
    print('* FP = ', fp)
    print('* FN = ', fn)
    print('*')
    print('************************')
    return tp, tn, fp, fn, recall, precision, num_predections


def get_k_mer_histogram(seq, k, printNonZero=False):
    seq = seq.upper()
    counter = Counter()
    bases = ['A', 'C', 'G', 'T']
    counter.update([seq[i:i + k] for i in range(len(seq) - k + 1)])

    # if printNonZero is True, output all kmers and thier frequency in seq:
    if printNonZero:
        for j in sorted(counter.keys()):
            print('{}\t{}'.format(j, counter[j]))

    kmer_and_values = [(''.join(p), 0) for p in itertools.product(bases, repeat=k)]
    histogram = []
    for (kmer, value) in kmer_and_values:
        if kmer in counter:
            histogram.append(counter[kmer])
        else:
            histogram.append(value)
    return histogram


def train_svm(data_path, k):
    with open(data_path, 'rb') as myf:
        data = pickle.load(myf)
    # Converting sequence data to histograms of k-mers
    print('Converting sequence data to histograms of k-mers...')
    X = []
    Y = []
    for (seq, data_label, _) in data:
        X.append(get_k_mer_histogram(seq, k))
        Y.append(data_label)

    size_of_test_set = int(len(Y) * TESTING_RATIO)

    clf = svm.SVC(gamma='scale')
    # print(' Started SVM Training...')
    clf.fit(X[:-size_of_test_set], Y[:-size_of_test_set])
    results = svm_testing_loss(clf, X[-size_of_test_set:], Y[-size_of_test_set:])
    return results


def compare_svm_on_k_list(data_path, k_list=range(1,7)):
    results = []
    for k in k_list:
        results.append(train_svm(data_path, k))
    return results


if __name__ == '__main__':
    # train_svm('data/all_data', K_MER_LEN)
    print(compare_svm_on_k_list('data/all_data'))

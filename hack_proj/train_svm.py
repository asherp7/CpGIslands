import itertools
from collections import Counter
import numpy as np
from sklearn import svm
import pickle
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.nan)

K_MER_LEN = 5
TESTING_RATIO = 0.1


# Results stats
def show_recall_precision_curve(recall, precision, title='SVM performance as function of K-mer length', show_grid=True, print_data=True):
    plt.figure()
    # lw = 2
    plt.plot(range(1, len(precision) + 1), precision, label='Precision')
    plt.plot(range(1, len(recall) + 1), recall, label='Recall')
    plt.xlabel('K')
    plt.ylabel('Percent')
    plt.legend()
    plt.title(title)
    # plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    # plt.xlim([0.0, 1.0])
    # plt.ylim([0.0, 1.05])
    # plt.xlabel('K')
    # # plt.ylabel('')
    # plt.title(title)
    # plt.legend(loc="lower right")
    # plt.grid(show_grid)
    #
    # # create the axis of thresholds (scores)
    # ax2 = plt.gca().twinx()
    # ax2.plot(fpr, thresholds, markeredgecolor='r', linestyle='dashed', color='g')
    # ax2.set_ylabel('Threshold', color='g')
    # ax2.set_ylim([thresholds[-1], thresholds[0]])
    # ax2.set_xlim([fpr[0], fpr[-1]])
    # if print_data:
    #     for i in range(fpr.size):
    #         print('FPR: %.3f, TPR: %.3f: Threshold = %.4f' % (fpr[i], tpr[i], thresholds[i]))
    plt.show()

def svm_testing_loss(clf, test_data, test_labels, k):
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
    print('* K-mer length = ', k)
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

    clf = svm.SVC(kernel='linear', gamma='scale')
    # print(' Started SVM Training...')
    clf.fit(X[:-size_of_test_set], Y[:-size_of_test_set])
    results = svm_testing_loss(clf, X[-size_of_test_set:], Y[-size_of_test_set:], k)
    bases = ['A', 'C', 'G', 'T']
    kmers = [''.join(p) for p in itertools.product(bases, repeat=k)]
    plot_coefficients(clf, kmers)
    return results


def compare_svm_on_k_list(data_path, k_list=range(1, 7)):
    results = []
    for k in k_list:
        results.append(train_svm(data_path, k))
    return results


def f_importances(coef, names):
    print(names)
    imp = coef
    imp, names = zip(*sorted(zip(imp, names)))
    plt.barh(range(len(names)), imp, align='center')
    plt.yticks(range(len(names)), names)
    plt.show()


def plot_coefficients(classifier, feature_names, top_features=15):
    coef = classifier.coef_.ravel()
    top_positive_coefficients = np.argsort(coef)[-top_features:]
    top_negative_coefficients = np.argsort(coef)[:top_features]
    top_coefficients = np.hstack([top_negative_coefficients, top_positive_coefficients])
    # create plot
    plt.figure(figsize=(15, 5))
    colors = ['red' if c < 0 else 'blue' for c in coef[top_coefficients]]
    plt.bar(np.arange(2 * top_features), coef[top_coefficients], color=colors)
    feature_names = np.array(feature_names)
    plt.xticks(np.arange(1, 1 + 2 * top_features), feature_names[top_coefficients], rotation=60, ha='right')
    plt.show()


if __name__ == '__main__':
    train_svm('data/all_data', K_MER_LEN)
    # results = compare_svm_on_k_list('data/data500/data_without_chr1')
    # recall = [x for (_,_,_,_,x,_,_) in results]
    # precision =[x for (_,_,_,_,_,x,_) in results]
    # recall.append(0.98185)
    # precision.append(0.99591)
    # show_recall_precision_curve(recall, precision)
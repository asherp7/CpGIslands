import re
import itertools
from utils import *

def list2dict(lst):
    d = {}
    for i in range(len(lst)):
        d[lst[i]] = i
    return d


class MarkovModel(object):
    UNKNOWN_VAL = -1

    def __init__(self, valid_states, order):
        self.order = order
        self.valid_states = valid_states
        self.valid_states_dict = list2dict(self.valid_states)
        self.states_num = len(self.valid_states)
        self.count_states = self._create_order_states(self.order + 1)
        self.count_states_dict = list2dict(self.count_states)
        self.base_states = self._create_order_states(self.order)
        self.base_states_dict = list2dict(self.base_states)

        self.log_prob_mat = np.log(np.ones(len(self.count_states), dtype=float) / len(self.count_states))
        self.log_transition_mat = np.log(np.ones((pow(self.states_num, self.order), self.states_num), dtype=float) / self.states_num)

    def get_prob_mat(self):
        return np.exp(self.log_prob_mat)

    def get_transition_mat(self):
        return np.exp(self.log_transition_mat)

    def _create_order_states(self, len):
        return list([''.join(p) for p in itertools.product(self.valid_states, repeat=len)])

    def _count_order_states(self, x, smooth=1):
        count_dict = {}
        for k in self.count_states:
            count_dict[k] = smooth

        for val in x:
            for i in range(len(val) - self.order):
                try:
                    count_dict[val[i:i + self.order + 1]] += 1
                except:
                    continue
        return count_dict

    def fit_transition(self, x):
        prob_mat = np.ones(len(self.count_states), dtype=float) / len(self.count_states)
        transition_mat = np.ones((pow(self.states_num, self.order), self.states_num), dtype=float) / self.states_num

        # Count all count states
        count_dict = self._count_order_states(x)

        # Count base states to normalize by
        norm_dict = {}
        for k in self.base_states:
            norm_dict[k] = 0
        for k, v in count_dict.items():
            norm_dict[k[:-1]] += v

        # Create probability mat for the count states
        for i, k in enumerate(self.count_states):
            prob_mat[i] = count_dict[k] / norm_dict[k[:-1]]

        # Create transition mat from base states
        for i, k in enumerate(self.base_states):
            for j, s in enumerate(self.valid_states):
                transition_mat[i, j] = count_dict[k + s] / norm_dict[k]

        self.log_prob_mat = np.log(prob_mat)
        self.log_transition_mat = np.log(transition_mat)

    def str_to_order_states(self, x):
        ret = []
        for val in x:
            curr_ret = []
            for i in range(len(val) - self.order):
                try:
                    curr_ret.append(self.count_states_dict[val[i:i + self.order + 1]])
                except:
                    curr_ret.append(MarkovModel.UNKNOWN_VAL)
            ret.append(np.array(curr_ret))
        return ret

    def get_ll(self, x):
        x_states = self.str_to_order_states(x)
        all_ll = []
        for arr in x_states:
            curr_ll = np.full_like(arr, fill_value=np.nan, dtype=float)
            for i in range(len(self.count_states)):
                curr_ll[arr == i] = self.log_prob_mat[i]
            self.nan_to_prev(curr_ll, np.mean(self.log_prob_mat))

            all_ll.append(curr_ll)
        return all_ll

    def nan_to_prev(self, arr, start_val):
        nan_idx = np.where(np.isnan(arr))[0]
        if nan_idx.size > 0:
            if np.isnan(arr[0]):
                arr[0] = start_val
            if nan_idx[0] == 0:
                if nan_idx.size == 1:
                    nan_idx = None
                else:
                    nan_idx = nan_idx[1:]
            if nan_idx is not None and nan_idx.size > 0:
                for i in nan_idx:
                    arr[i] = arr[i - 1]


def print_prob_mat(markov_model: MarkovModel):
    prob_mat = markov_model.get_prob_mat()
    for i, k in enumerate(markov_model.count_states):
        print('%s: %f' % (k, prob_mat[i]))


def print_transition_mat(markov_model: MarkovModel):
    states_str = ''
    transition_mat = markov_model.get_transition_mat()
    for j in range(transition_mat.shape[1]):
        states_str += '\t%s' % markov_model.valid_states[j]
    print(states_str)
    for i in range(transition_mat.shape[0]):
        print_str = markov_model.base_states[i]
        for j in range(transition_mat.shape[1]):
            print_str += '\t%.3f' % transition_mat[i,j]
        print(print_str)


if __name__ == "__main__":
    m = MarkovModel(['A', 'B'], 1)
    train = ['12AAG', 'AABBAA', 'ABABAAAAAABBBBBNNABABABA', 'AAXXABBBAA', 'ABBBBBBXBBBAABBBBBBBBBBAAAAAAAAAAAAAAAAAAAAAAAAA']
    print_prob_mat(m)
    print_transition_mat(m)
    m.fit_transition(train)
    print('-----------------------')
    print_prob_mat(m)
    print_transition_mat(m)

    ll = m.get_ll(train)
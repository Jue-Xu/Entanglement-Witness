# Various utility functions used by the scripts

import itertools
import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
import math
from cmath import cos, sin, exp, pi, sqrt
import random

from qutip import *

from sklearn import svm, datasets
from sklearn.inspection import DecisionBoundaryDisplay
from sklearn.svm import SVC
from sklearn.feature_selection import RFE

########################################
""" PUBLIC FUNCTIONS USED IN NOTEBOOK """
########################################
color_config = {'entangled': 'c', 'separable': 'm', 'unfaithful': 'y', 'ghz': 'r', 'w': 'g', 'cluster': 'b', 'bell': 'k'}
color_dict = {0: 'c', 1: 'm', 2: 'y'}

pauli_operators = [qeye(2), sigmax(),sigmay(),sigmaz()]
pauli_str = ['I', 'X', 'Y', 'Z']

""" Random states creation """
def generate_rand_product_state(n, m, noise_limit=0, is_pure=True):
    random_white_noise_p = [random.random() * noise_limit for i in range(m)]
    dim = [2 for i in range(n)]

    if is_pure:
        print('generate_rand_pure_product_state+noise:', n, 'qubits,', m,
              'samples')
        return [
            ket2dm(tensor([rand_ket(2) for j in range(n)])) * (1 - p) +
            p / 2**n * identity(dims=dim) for p in random_white_noise_p
        ]
    else:
        print('generate_rand_product_density_matrix:', n, 'qubits,', m,
              'samples')
        # print(rand_dm(2))
        # print(tensor([rand_dm(2),rand_dm(2),rand_dm(2)]))
        return [tensor([rand_dm(2) for j in range(n)]) for i in range(m)]


def generate_rand_product_density(n, m, noise_limit):
    print('generate_rand_product_density:', n, 'qubits,', m, 'samples')
    # print(rand_dm(2))
    # print(tensor([rand_dm(2),rand_dm(2),rand_dm(2)]))
    return [tensor([rand_dm(2) for j in range(n)]) for i in range(m)]


# generate_rand_product_density(3,1,0)


def generate_noisy_biseparable(n_A, n_B, m, noise_limit):
    random_white_noise_p = [random.random() * noise_limit for i in range(m)]
    # print(tensor(rand_dm([2]),rand_dm([2,2])))
    dim_AB = [2 for i in range(n_A + n_B)]
    dim_A = [[2 for i in range(n_A)], [2 for i in range(n_A)]]
    dim_B = [[2 for i in range(n_B)], [2 for i in range(n_B)]]
    # print(rand_dm(2))
    # print(identity([2 for i in range(3)]))
    # print(rand_dm(N=2,dims=[[2],[2]]))
    # print(rand_dm(N=4,dims=[[2,2],[2,2]]))

    return [
        tensor(rand_dm(N=2**n_A, dims=dim_A), rand_dm(N=2**n_B, dims=dim_B)) *
        (1 - p) + p / 2**(n_A + n_B) * identity(dim_AB)
        for p in random_white_noise_p
    ]


# generate_noisy_biseparable(1,2,10,1/3)


def generate_bell_noisy_density(m, kind, noise_limit):
    # p_limit = 1/3
    # singlet_state()

    random_white_noise_p = [random.random() * noise_limit for i in range(m)]

    return [
        ket2dm(bell_state(state=kind)) * (1 - p) + p / 4 * identity([2, 2])
        for p in random_white_noise_p
    ]


# def generate_two_qubit_entangled_pure_state(m):


def generate_bell_like_pure_state(m):
    # permute(order)
    theta_list = [random.random() * 2 * pi for i in range(m)]
    phi_list = [random.random() * pi for i in range(m)]

    a_list = [cos(theta) for theta in theta_list]
    b_list = [
        sin(theta) * exp(phi * 1j) for theta in theta_list for phi in phi_list
    ]
    # x = a* basis(4, 0) + b* basis(4, 3)
    # print(x.norm())
    # print(x)
    # print(x.isket)

    return [
        ket2dm(a * tensor(basis(2, 0), basis(2, 0)) +
               b * tensor(basis(2, 1), basis(2, 1))) for a in a_list
        for b in b_list
    ]


def generate_noisy_ghz_ensemble(n, m, noise_limit):
    print('generate_noisy_ghz_ensemble:', n, 'qubits,', m, 'samples')
    return [
        ket2dm(ghz_state(N=n)) * (1 - p_noise) + p_noise /
        (2**n) * qeye([2 for j in range(n)])
        for p_noise in [random.random() * noise_limit for i in range(m)]
    ]


def generate_noisy_w_ensemble(n, m, noise_limit):
    return [
        ket2dm(w_state(N=n)) * (1 - p_noise) + p_noise /
        (2**n) * qeye([2 for j in range(n)])
        for p_noise in [random.random() * noise_limit for i in range(m)]
    ]


# bell_like_pure_state = generate_bell_like_pure_state(10, False)


def const_label(c, m):
    return [c for i in range(m)]


# generate_rand_product_state(4,3)


def ppt_criterion(rho):
    # Necessary and sufficient condition (2-qubit): Positive Partial Transpose (PPT)
    if len(rho.dims[0]) == 2:
        # print(rho.dims[0])
        rho_out = partial_transpose(rho, [0, 1])
    elif len(rho.dims[0]) == 3:
        # print(rho.dims[0])
        rho_out = partial_transpose(rho, [1, 1, 0])
    else:
        return len(rho.dims[0])
    smallest_eigenval = rho_out.eigenenergies(sort='low', eigvals=1)
    # print(rho_out.eigenenergies(sort='low',eigvals=1))
    # print(smallest_eigenval)
    return smallest_eigenval


def generate_two_qubit_random_state_PPT(m, plot=False):
    rand_dm_2 = [rand_dm(N=4, dims=[[2, 2], [2, 2]]) for i in range(m)]
    # print(ppt_criterion(random.choice(rand_dm_2)))
    entangled = [rho for rho in rand_dm_2 if ppt_criterion(rho) < 0]
    # entangled = [filter(lambda p: ppt_criterion(p) < 0, rand_dm_2)]
    # print(entangled)

    # for product pure states, the smallest eigenvalue is a very very small negative one
    separable = [rho for rho in rand_dm_2 if ppt_criterion(rho) >= 0]
    print('# entangled state:', len(entangled), '; # separable state:',
          len(separable))
    # print()

    # smallest_eigen_list_random = np.array([ppt_criterion(state) for state in rand_dm_2]).flatten()

    if plot:
        fig, ax = plt.subplots(figsize=(6, 4))

        # print(np.array(smallest_eigen_list_entangled).flatten())
        # ax.hist(smallest_eigen_list_random)
        ax.hist(np.array([ppt_criterion(state)
                          for state in entangled]).flatten(),
                color=color_dict[0])
        ax.hist(np.array([ppt_criterion(state)
                          for state in separable]).flatten(),
                color=color_dict[1])
        # ax.text(4, 0.4, r'$$')
        ax.set_ylabel('Samples')
        ax.set_xlabel('smallest eigenvalue of partial transpose')
        ax.set_title('2-qubit random density matrix (test PPT criterion)')

    # plt.savefig('two_qubit_PPT_hist.png', dpi=400)

    return [entangled, separable]

#####################################################
# machine learning
#####################################################
def construct_training_dataset(states_labels,operators,verbose=False):
    print("------- construct_training_dataset -------")
    states_labels = np.array(states_labels,dtype=object)
    all_states = states_labels[:,0].flatten()
    # print(all_states)
    y = states_labels[:,1].flatten()
    # print(y)
    # y = sum(states_labels[:,1],[])
    n_sample = len(all_states)

    # evaluate features - expectation values
    features = np.array([ expect(operators, state) for state in all_states ])

    np.random.seed(1)
    order = np.random.permutation(n_sample)
    # print(order)
    features = features[order]
    y = y[order].astype(int)
    if verbose:
        print(y)
    # print(features)

    print('number of samples:', n_sample, '; number of labels:', len(y), '; dimension:', shape(features))
    print("------------------- end -------------------")
    return (features, y)


def plot_ranking(x_labels, ranking):
    fig_rank, ax_rank = plt.subplots(figsize=(6, 4))
    barlist = ax_rank.bar(x_labels, ranking, edgecolor="k")
    for index, item in enumerate(ranking):
        if item == max(ranking):
            barlist[index].set_color('gray')
        if item == min(ranking):
            barlist[index].set_color('red')

    ax_rank.set_ylabel('ranking')
    ax_rank.set_xlabel('Pauli operator')
    # ax_rank.legend(['','',''])
    plt.savefig('feature_rank.png', dpi=300)


def plot_feature_space(X, y, filter):
    fig_feature, ax_feature = plt.subplots(figsize=(6, 4))
    # ax_feature.scatter(X_train[:, filter][:, 0], X_train[:, filter][:, 1], c=y_train, zorder=10, cmap=plt.cm.Paired, edgecolor="k", s=20)
    # https://stackoverflow.com/questions/47006268/matplotlib-scatter-plot-with-color-label-and-legend-specified-by-c-option

    # https://matplotlib.org/stable/gallery/color/named_colors.html
    for g in np.unique(y):
        ix = np.where(y == g)
        ax_feature.scatter(X[:, filter][ix, 2],
                           X[:, filter][ix, 1],
                           c=color_dict[g],
                           edgecolor="k",
                           s=30)
    ax_feature.set_ylabel('Select feature #2')
    ax_feature.set_xlabel('Select feature #1')
    ax_feature.legend(['entangled', 'separable'], loc='upper right')
    plt.savefig('feature_space.png', dpi=300)


def my_svm(X, y, size_test, kernel, legend, rfe=False, to_features=3):
    ################# SVM training ####################
    print(
        "=========================== SVM summary start ============================"
    )
    print("size of training set:", len(X), "; size of testing set:", size_test)
    # print("size of testing set:", size_test)
    n_sample = len(X)
    X_train = X[:int(0.9 * n_sample)]
    # print(X_train)
    y_train = y[:int(0.9 * n_sample)]
    X_test = X[int(0.9 * n_sample):]
    y_test = y[int(0.9 * n_sample):]

    # we create an instance of SVM and fit out data.
    if kernel == 'linear':
        # linear kernel
        print("kernel method: linear kernel")
        if rfe:
            print("---- recursive feature elimination ----")
            estimator = SVC(kernel="linear", C=1)
            clf = RFE(estimator=estimator,
                      n_features_to_select=to_features,
                      step=1)
            clf.fit(X_train, y_train)
            filter = clf.support_
            # print('feature filter:', filter)
            ranking = clf.ranking_.reshape(X_train[0].shape)
            print('feature ranking:', ranking)
            print('----------------------------------------')
            # print(y_train)
            plot_ranking(two_pauli_label, ranking)
            plot_feature_space(X_train, y_train, filter)
        else:
            clf = svm.SVC(kernel='linear')
            # clf = svm.SVC(kernel=my_kernel)
            # # clf.get_params()
            clf.fit(X_train, y_train)
            # print('score:', clf.score)
            print('coef0:', clf.coef0)
            print('coef_:', clf.coef_)
            print('intercept_:', clf.intercept_)
    else:
        # kernel
        print("kernel method: rbf")
        clf = svm.SVC()
        clf.fit(X_train, y_train)

    train_score = clf.score(X_train, y_train)
    test_score = clf.score(X_test, y_test)
    print('score (train):', train_score)
    print('score (test):', test_score)
    # print(clf.score)

    fig, ax = plt.subplots(figsize=(6, 4))

    ################# test/prediction ####################
    prediction_test = clf.predict(X_test)
    decision_test = clf.decision_function(X_test)
    # # print(prediction)
    # print("accuracy of prediction 0 (bell entangled):", sum((1-prediction_0))/len(prediction_0))
    ax.hist(decision_test, bins=25)
    ax.axvline(0.0, ls="--", color="r")

    # test_2 = generate_bell_like_pure_state(10)
    test_2 = generate_two_qubit_random_state_PPT(size_test)[1]
    # test_2 = generate_rand_product_state(2,size_test)
    # test_2 = generate_bell_noisy_density(size_test,'10',1/3)
    feature_2 = [expect(two_pauli, state) for state in test_2]
    prediction_2 = clf.predict(feature_2)
    # print(prediction_2)
    decision_2 = clf.decision_function(feature_2)
    ax.hist(decision_2, alpha=0.7, bins=25)
    print("accuracy of prediction (other entangled):",
          sum(prediction_2) / len(prediction_2))

    # ax.hist([expect(state, bell_inequality) for state in test_0 ])
    # ax.hist([expect(state, bell_inequality) for state in [ ket2dm(ket) for ket in test_1 ] ])
    # ax.hist([expect(state, bell_inequality) for state in test_2 ])

    print(
        "============================= SVM summary end =============================="
    )
    ax.set_ylabel('Number of samples')
    ax.set_xlabel('Expectation value')
    ax.set_title('2-qubit')
    ax.legend(legend, loc='upper right')
    plt.savefig('two_qubit_hist.png', dpi=300)

    # return clf
    return (train_score, test_score)



def plot_score(size_list, train_score_list, test_score_list):
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(size_list, train_score_list,"o-")
    ax.plot(size_list, test_score_list,"o-")
    ax.legend(['train','test'], loc='upper right')
    ax.set_ylabel('Score (accuray)')
    ax.set_xlabel('Number of samples')
    plt.savefig('two_qubit_scores.png', dpi=300)
    

def plot_expectation_hist(ax, expectation_list, legends, title=''):
    ax.axvline(0.0, ls="--", color="gray")
    for expectation in expectation_list:
        ax.hist(expectation, alpha=0.7, edgecolor="k",linewidth=0.5)
    ax.legend(legends, loc='upper left', prop ={'size': 9})
    ax.set_ylabel('Number of samples')
    ax.set_xlabel('Expectations of different entanglement witness')
    ax.text(0.04, 0.62, title, transform=ax.transAxes)
    # ax.set_title('3-qubit case '+title)
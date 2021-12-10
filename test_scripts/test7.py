# -*- coding: utf-8 -*-
# Created on Mon Nov  6 21:39:46 2017"

__author__ = "Siqi Miao"
__contributor__ = "Mohit M."

import numpy as np
import random
from numpy.linalg import norm, cond

eps = 1.0e-10


def run(func):
    """
    test7.py
    Feasible initialization test.
    indices of iB, iN start with 1

    :param func: Construction of simplex_init function
    :return: Test Success/Fail
    """
    test_name = "TEST-7"
    assert callable(func), "Argument must be a callable function"
    assert func.__name__ == "simplex_init", test_name + " checks the simplex_init function. Incorrect callable passed!"

    print("************************************")
    print("\t{}".format(test_name))
    print("************************************")
    test_success = True

    # first form an invertible matrix
    R = np.matrix([[4, 1, 1],
                   [-1, 2, 1],
                   [1, 1, -1]], dtype=np.float64)

    # form a vector b which is in the span of R
    b = R * abs(randn(3, 1))

    B = np.matrix([[1, 1, 1],
                   [1, 1, 0],
                   [1, 0, 0]], dtype=np.float64)
    A = np.hstack((R, B))

    # form a random permutation
    p = list(range(0, 6))
    random.shuffle(p)
    A = A[:, p]

    # c doesn't matter for this test
    c = np.matrix([[-1, -1, -1, -1, -1, -1]], dtype=np.float64)
    [istatus, iB, iN, xB] = func(A, b, c)

    X = np.zeros((6, 1), dtype=np.float64)
    X[[(b1 - 1) for b1 in iB]] = xB

    if istatus != 0:
        print('>> istatus wrong!')
        test_success = False

    # test feasibility
    if norm(A * X - b) > eps:
        print('>> NOT FEASIBLE!!!')
        test_success = False

    if min(X) < 0:
        print('>> NOT FEASIBLE!!!')
        test_success = False

    # test that we have a basis
    if (1 / cond(A[:, [(b - 1) for b in iB]])) > 1.0e6:
        print('>> NOT BASIC!!!')
        test_success = False

    return test_success

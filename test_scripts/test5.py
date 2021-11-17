# -*- coding: utf-8 -*-
# Created on Mon Nov  6 21:39:46 2017"

__author__ = "Siqi Miao"
__contributor__ = "Mohit M."

import numpy as np
from numpy.linalg import norm


def run(func):
    """
    test5.py
    Unboundedness test.
    indices of iB, iN start with 1

    :param func: Construction of simplex_step function
    :return: Test Success/Fail
    """
    test_name = "TEST-5"
    assert callable(func), "Argument must be a callable function"
    assert func.__name__ == "simplex_step", test_name + " checks the simplex_step function. Incorrect callable passed!"

    print("************************************")
    print("\t{}".format(test_name))
    print("************************************")
    test_success = True

    # start with a tableau form
    A1 = np.matrix([[-1, 1, 2],
                    [-1, 1, 1],
                    [0, 1, 1]], dtype=np.float64)

    A = np.hstack((np.eye(3), A1))

    b = np.matrix([[1],
                   [2],
                   [3]], dtype=np.float64)

    iB = [1, 2, 3]
    iN = [4, 5, 6]
    xB = np.matrix(np.copy(b))
    c = np.matrix([[0, 0, 0, -1, 2, 1]], dtype=np.float64)

    # form an invertible matrix B and modify the problem
    B = np.matrix([[4, 1, 0],
                   [1, -2, -1],
                   [1, 2, 4]], dtype=np.float64)
    A = B * A
    b = B * b

    # modify c
    N = A[:, [index_N - 1 for index_N in iN]]
    c1 = np.matrix([[1, 1, 0]], dtype=np.float64)
    c2 = c[:, (4 - 1):6] + c1 * B.I * N
    c = np.hstack((c1, c2))

    # take a step.
    irule = 0
    [istatus, iB, iN, xB] = func(A, b, c, iB, iN, xB, irule)

    X = np.zeros((6, 1), dtype=np.float64)
    X[[(b - 1) for b in iB]] = xB

    if (istatus != 16):
        print('>> INCORRECT ISTATUS!')
        test_success = False

    return test_success

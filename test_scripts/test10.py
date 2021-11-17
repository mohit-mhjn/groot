# -*- coding: utf-8 -*-
# Created on Mon Nov  6 21:39:46 2017"

__author__ = "Siqi Miao"
__contributor__ = "Mohit M."

import numpy as np
from .simplex_test import simplex_test

eps = 1.0e-10


def run(func):
    """
    test10.py
    An infeasible input to simplex method.
    indices of iB, iN start with 1

    :param func: Construction of simplex_method function
    :return: Test Success/Fail
    """
    test_name = "TEST-10"
    assert callable(func), "Argument must be a callable function"
    assert func.__name__ == "simplex_method", test_name + " checks the simplex_method function. Incorrect callable passed!"

    print("************************************")
    print("\t{}".format(test_name))
    print("************************************")
    test_success = True

    # first form an invertible matrix
    R = np.matrix([[4, 1, 1],
                   [1, 2, 1],
                   [1, 1, 1]], dtype=np.float64)

    # form a vector b which is in the span of R
    b = R * np.matrix([[1],
                       [-4],
                       [-1]], dtype=np.float64)

    B = np.matrix([[1, 1, 1],
                   [1, 1, 0],
                   [1, 0, 0]], dtype=np.float64)

    A = np.hstack((R, B))
    c = np.matrix([[-2, 1, 1, -1, -1, -1]], dtype=np.float64)

    irule = 1
    [istatus, X, eta, iB, iN, xB] = func(A, b, c, irule)

    if istatus != 4:
        print('>> istatus is wrong')
        test_success = False

    return test_success

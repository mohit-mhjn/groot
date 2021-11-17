# -*- coding: utf-8 -*-
# Created on Mon Nov  6 21:39:46 2017"

__author__ = "Siqi Miao"
__contributor__ = "Mohit M."

import numpy as np
from numpy.linalg import norm


def run(func):
    """
    test6.py
    First initialization test: infeasible problem.
    indices of iB, iN start with 1

    :param func: Construction of simplex_init function
    :return: Test Success/Fail
    """
    test_name = "TEST-6"
    assert callable(func), "Argument must be a callable function"
    assert func.__name__ == "simplex_init", test_name + " checks the simplex_init function. Incorrect callable passed!"

    print("************************************")
    print("\t{}".format(test_name))
    print("************************************")
    test_success = True

    A = np.matrix([[1, 1, 1, 2, 1, 3],
                   [1, 1, 0, 2, 2, 2],
                   [1, 0, 0, 12, 1, 1]], dtype=np.float64)

    b = np.matrix([[-1],
                   [3],
                   [-1]], dtype=np.float64)

    c = np.matrix([[-1, -1, -1, -1, -1, -1]], dtype=np.float64)

    [istatus, iB, iN, xB] = func(A, b, c)

    if (istatus != 16):
        print('>> istauts WRONG!!!!')
        test_success = False

    A = np.matrix(np.copy(-A))
    c = np.matrix(np.copy(c))

    [istatus, iB, iN, xB] = func(A, b, c)

    if (istatus != 16):
        print('>> istauts WRONG!!!!')
        test_success = False

    return test_success

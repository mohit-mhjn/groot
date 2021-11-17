# -*- coding: utf-8 -*-
# Created on Mon Nov  6 21:39:46 2017"

__author__ = "Siqi Miao"
__contributor__ = "Mohit M."

import numpy as np
from numpy.linalg import norm


def run(func):
    """
    test2.py

    Test that the simplex method code takes a correct step
    when the form is slightly more complicated.

    indices of iB, iN start with 1
    :param func: Construction of simplex_step function
    :return: Test Success/Fail
    """
    test_name = "TEST-2"
    assert callable(func), "Argument must be a callable function"
    assert func.__name__ == "simplex_step", test_name + " checks the simplex_step function. Incorrect callable passed!"

    print("************************************")
    print("\t{}".format(test_name))
    print("************************************")
    test_success = True

    A = np.matrix([[-4, 1, 0, -3, -3, -5],
                   [1, -2, -1, -2, -2, 3],
                   [1, 2, 4, 7, 7, -1]], dtype=np.float64)

    b = np.matrix([[-2],
                   [-6],
                   [17]], dtype=np.float64)

    iB = [1, 2, 3]
    iN = [4, 5, 6]
    xB = np.matrix([[1], [2], [3]], dtype=np.float64)
    c = np.matrix([[1, 1, 1, 2, 5, 1]], dtype=np.float64)

    # take a step
    irule = 0
    [istatus, iB, iN, xB] = func(A, b, c, iB, iN, xB, irule)

    X = np.zeros((6, 1), dtype=np.float64)
    X[[(b - 1) for b in iB]] = xB

    if istatus != 0:
        print('>> INCORRECT ISTATUS!')
        test_success = False

    if norm(X - np.matrix([[0], [1], [2], [1], [0], [0]])) > 1e-10:
        print('>> INCORRECT STEP!')
        test_success = False

    if norm(np.array(sorted(iN)) - np.array([1, 5, 6])) > 1e-10:
        print('>> iN incorrect!')
        test_success = False

    if norm(np.array(sorted(iB)) - np.array([2, 3, 4])) > 1e-10:
        print('>> iB incorrect!')
        test_success = False

    return test_success

# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 16:04:29 2017

@author: Siqi Miao
"""

import numpy as np
from numpy.linalg import norm
from functools import wraps
from pprint import pprint


def simplex_test(A, b, c, iB, xB):
    """
    Test the feasibility and optimality of a basic vector for the linear program:

        min:          cx
        subject to:   Ax=b                                    (1)
                      x>= 0

    where A is an (m,n) matrix of rank m.

    Input Parameters:

    A - (m,n) constraint matrix
    b - (m,1) vector appearing in the constraint equation above
    c - (1,n) vector giving the coefficients of the objective function

    iB - (1,m) integer vector specifying the indices of the basic
         variables at the beginning of the simplex step
    xB - (m,1) vector specifying the values of the basic
         variables at the beginning of the simplex step

    Output Parameters:

    X - vector of length n which contains both the basic and
       nonbasic values
    eta - value of the objective function at X
    isfeasible - integer parameter indicating whether the given basic vector is feasible or not;
         isfeasible = 0  means the vector is infeasible
         isfeasible = 1  means the vector is feasible
    zN - the value of the dual basic variables corresponding to the given basic feasible vector
    indices of iB, iN start with 1
    """
    iB = [i - 1 for i in iB]

    [m, n] = A.shape
    eps = 1.0e-12

    X = []
    eta = []
    isfeasible = []
    isoptimal = []

    X = np.zeros((n, 1))
    X[iB] = xB
    eta = c * X

    err = norm(A * X - b)

    isfeasible = 0
    if (err < eps) and min(X) >= -eps:
        isfeasible = 1

    temp = list(range(n))
    iN = []
    for each in temp:
        if each not in iB:
            iN.append(each)

    Cb = c[:, iB]
    Cn = c[:, iN]
    B = A[:, iB]
    N = A[:, iN]

    Binv = B.I

    ctilde = []
    for i in range(len(iN)):
        ctilde.append(Cn[:, i] - Cb * Binv * A[:, iN[i]])

    if not (min(ctilde) >= -eps):
        isoptimal = 0

    else:
        isoptimal = 1

    zN = ctilde

    return [X, eta, isfeasible, isoptimal, zN]


def protect_step(func):
    @wraps(func)
    def argument_validator(*args, **kwargs):
        # simplex_step(A,b,c,iB,iN,xB,irule)
        # assert isinstance(A, np.matrix)
        [istatus, iB, iN, xB] = func(*args, **kwargs)
        assert isinstance(istatus, int), "return: istatus must be int type"
        assert isinstance(iB, list), "return: iB must be list type"
        assert isinstance(iN, list), "return: iN must be list type"
        assert isinstance(xB, np.matrix), "return: xB must be np.matrix type"
        return [istatus, iB, iN, xB]

    return argument_validator


def protect_init(func):
    @wraps(func)
    def argument_validator(*args, **kwargs):
        # simplex_step(A,b,c,iB,iN,xB,irule)
        # assert isinstance(A, np.matrix)
        [istatus, iB, iN, xB] = func(*args, **kwargs)
        assert isinstance(istatus, int), "return: istatus must be int type"
        assert isinstance(iB, list), "return: iB must be list type"
        assert isinstance(iN, list), "return: iN must be list type"
        assert isinstance(xB, np.matrix), "return: xB must be np.matrix type"
        return [istatus, iB, iN, xB]

    return argument_validator


def protect_method(func):
    @wraps(func)
    def argument_validator(*args, **kwargs):
        # simplex_step(A,b,c,iB,iN,xB,irule)
        # assert isinstance(A, np.matrix)
        [istatus, X, eta, iB, iN, xB] = func(*args, **kwargs)
        assert isinstance(istatus, int), "return: istatus must be int type"
        assert isinstance(iB, list), "return: iB must be list type"
        assert isinstance(iN, list), "return: iN must be list type"
        assert isinstance(xB, np.matrix), "return: xB must be np.matrix type"
        return [istatus, X, eta, iB, iN, xB]

    return argument_validator


# todo: better detailing for obtained and expected -
#   derive formatting and comparison based on input instance type
class groot(object):
    def __init__(self, test_name, func_name):
        """
        Multipurpose Debugger and Tester for the use cases
        Mohit - still working on this passively
        :param test_name:
        :param func_name:
        """
        self.test_name = test_name
        self.func_name = func_name
        self.super_debug = False

    def match_iB(self):
        pass

    def match_iN(self):
        pass

    def match_istatus(self):
        pass

    def check_positive(self):
        pass

    def check_feasible(self, expression):
        pass

    def optimal(self):
        pass

    def pretty_print_comparator(self):
        pass


def pretty_print_debug(expected, obtained, test_name, func_name, item, fail=True):
    """
    Mohit - still working on this passively
    :param expected:
    :param obtained:
    :param test_name:
    :param func_name:
    :param item:
    :param fail:
    :return:
    """
    print("-" * 40)
    if fail:
        print("Mismatch found in: {} for function: {} and return: {}".format(test_name, func_name, item))
    else:
        print("Match found in: {} for function: {} and return: {}".format(test_name, func_name, item))
    print("-" * 40)
    print("*" * 10, "EXPECTED VALUE >>>>")
    pprint(expected)
    print("*" * 10, "OBTAINED VALUE >>>>")
    pprint(obtained)
    print("-" * 40)
    return None


def check_failure(obtained, expected, op, test_name, func, func_return, debug):
    super_debug = False
    if op(obtained, expected):
        print('>> INCORRECT return value for : {}'.format(func_return))
        if debug:
            pretty_print_debug(expected, obtained, test_name, func, func_return)
        return False
    elif super_debug:
        pretty_print_debug(expected, obtained, test_name, func, func_return, fail=False)

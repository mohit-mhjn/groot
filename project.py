from test_scripts.simplex_test import protect_step, protect_init, protect_method
import numpy as np
from numpy.linalg import inv


@protect_step
def simplex_step(A, b, c, iB, iN, xB, Binv, irule):
    """
    Take a single simplex method step for the linear program

    min:    c*x
    ST:     Ax=b
            x>=0,

    where A is an (m,n) matrix.

    That is, given a basic feasible vector vector described by the variables
    iB,iN,XB,

    Input Parameters:

    A - (n,m) constraint matrix
    b - (m,1) POSITIVE vector appearing in the constraint equation above
    c - (1,n) vector giving the coefficients of the objective function

    iB - (1,m) integer vector specifying the indices of the basic variables at the beginning of the simplex step
    iN - (1,n-m) integer vector specifying the indices of the nonbasic variables at the beginning of the simplex step
    xB - (m,1) vector specifying the values of the basic variables at the beginning of the simplex step

    irule - integer parameter specifying which pivot rule to use:
        irule = 0 indicates that the smallest coefficient rule should be used
        irule = 1 indicates that Bland's rule should be used

    Output Parameters:

    istatus - integer parameter reporting the condition of the
       istatus = 0  indicates normal simplex method step completed
       istatus = 16 indicates the program is unbounded
       istatus = -1 indicates an optimal feasible vector has been found

    iB - integer vector specifying the m indices of the basic variables after the simplex step
    iN - integer vector specifying the n-m indices of the nonbasic variables after the simplex step
    xB - vector of length m specifying the values of the basic variables after the simplex step
    """
    istatus = 0

    # @@ Write program here to update values for [istatus, iB, iN, xB] >>
    # [Cite: BT Linear Optimization : Revised Simplex - Pag 114]

    # STEP 1. Compute the reduced costs using following >>
    # zj = cB^T * B^-1 * A.j
    # cj - zj

    cB = np.matrix([[c[0, i - 1]] for i in iB])
    entering_variable_candidates = {}

    for j in iN:
        cj_minus_zj = (c[0, j - 1] - cB.transpose() * Binv * A[:, j - 1])[0, 0]
        if cj_minus_zj < 0:
            entering_variable_candidates[j] = cj_minus_zj

    # Step 1.1 >> Check if everything was optimal/non-negative reduced costs >>
    if not entering_variable_candidates:
        istatus = -1
        return [istatus, iB, iN, xB, Binv]

    # There is something negative => There is at least one entering variable >>
    # STEP 2 :: Chose the leaving variable rule:
    if irule == 0:
        # smallest coefficient rule >> Dantzig's rule (choose most negative reduced costs)
        entering_var_index = min(entering_variable_candidates, key=entering_variable_candidates.get)
    elif irule == 1:
        # Bland's Rule - function picks smallest index for entering variable
        entering_var_index = min(i for i in entering_variable_candidates.keys())
    else:
        raise AssertionError("Invalid irule!")

    # STEP 3: Do a minimum ratio test >>
    # RHS from tableau >>
    b_hat = Binv * b
    b_inv_aj = Binv * A[:, entering_var_index - 1]

    ratio = 99999
    leaving_var_index = -1
    leaving_var_position = -1

    for position, b_var_indx in enumerate(iB):
        if b_inv_aj[position] > 0:
            temp_ratio = b_hat[position, 0] / b_inv_aj[position, 0]
            if temp_ratio < ratio:
                ratio = temp_ratio
                leaving_var_index = b_var_indx
                leaving_var_position = position

    # If no ratio was positive then >> Minimum ratio test failed and solution is unbounded >>
    if leaving_var_index == -1:
        istatus = 16
        return [istatus, iB, iN, xB, Binv]

    # If not unbounded >> then use the leaving variable index and do the pivot
    # STEP 4: Decided the leaving and entering variable
    xB = np.matrix([[xB[position, 0] - ratio * b_inv_aj[position, 0]] for position, _ in enumerate(iB)])
    xB[leaving_var_position, 0] = ratio

    iB[iB.index(leaving_var_index)] = entering_var_index
    iN[iN.index(entering_var_index)] = leaving_var_index

    # STEP 5: Find the updated B inverse -- revised simplex way
    Binv_concat_y = np.hstack([Binv, b_inv_aj])
    updated_B_inv_concat_y = np.matrix(np.zeros(Binv_concat_y.shape))
    for col in range(Binv_concat_y.shape[1]):
        updated_B_inv_concat_y[leaving_var_position, col] = Binv_concat_y[leaving_var_position, col] / Binv_concat_y[
            leaving_var_position, -1]
    for row in range(Binv_concat_y.shape[0]):
        if not row == leaving_var_position:
            for col in range(Binv_concat_y.shape[1]):
                updated_B_inv_concat_y[row, col] = Binv_concat_y[row, col] - (
                        Binv_concat_y[row, -1] / updated_B_inv_concat_y[
                    leaving_var_position, -1]) * updated_B_inv_concat_y[leaving_var_position, col]

    Binv = np.delete(updated_B_inv_concat_y, -1, axis=1)
    istatus = 0
    return [istatus, iB, iN, xB, Binv]


@protect_init
def simplex_init(A, b, c, irule=1):
    """
    Attempt to find a basic feasible vector for the linear program

    max:    c*x
    ST:     Ax=b
            x>=0,

    where A is a (m,n) matrix.

    Input Parameters:

    A - (n,m) constraint matrix
    b - (m,1) vector appearing in the constraint equation above
    c - (1,n) vector giving the coefficients of the objective function

    Output Parameters:

    istatus - integer parameter reporting the condition of the
       istatus = 0  indicates a basic feasible vector was found
       istatus = 4  indicates that the initialization procedure failed
       istatus = 16  indicates that the problem is infeasible

    iB - integer vector of length m specifying the indices of the basic variables
    iN - integer vector of length n-m specifying the indices of the nonbasic variables
    xB - vector of length m specifying the values of the basic variables
    """
    istatus = None
    iB = None
    iN = None
    xB = None

    # @@ Write program here to assign appropriate values to [istatus, iB, iN, xB] >>>
    rows = A.shape[0]
    columns = A.shape[1]

    # Check if any of the b's were negative >>
    # If that happens, then multiply the constraint by -1
    for i in range(rows):
        if b[i, 0] < 0:
            for j in range(columns):
                A[i, j] *= -1
            b[i, 0] *= -1

    # PREPARE ARTIFICIAL VARIABLES >>
    artificial_vars_matrix = np.eye(rows)
    artificalA = np.hstack([A, artificial_vars_matrix])

    # c = np.hstack((c, len(A[:,0])))
    # cost vector for big_M row >>
    costM = np.matrix([[0 for _ in range(columns)] + [1 for _ in range(rows)]], dtype=float)

    # Row Operations to make all artificial variables basic for initialization
    # for i in range(rows):
    #     costM -= artificalA[i, :]

    iB = [columns + i + 1 for i in range(rows)]
    iN = [i + 1 for i in range(columns)]
    iA = list(iB)  # artificial variable indices

    # initializing B-inverse and xB
    Binv = np.eye(rows)
    xB = b
    istatus = 4
    termination_istatus = set([-1, 16])
    max_iterations = 99  # >>> 300 times of number of columns
    counter = 0

    # iterating to optimality
    while (istatus not in termination_istatus and counter < max_iterations):
        # calling simplex_step for variable values after one iteration
        [istatus, iB, iN, xB, Binv] = simplex_step(artificalA, b, costM, iB, iN, xB, Binv, irule=irule)
        counter += 1

    if counter >= max_iterations or istatus == 16:
        # function terminates here with unknown occurrence
        return [4, iB, iN, xB]

    artifical_vars_in_basis = set(iA).intersection(set(iB))

    # Test Level 1:
    if not artifical_vars_in_basis:
        iN = [i for i in iN if i not in set(iA)]
        return [0, iB, iN, xB]
    else:
        for some_var in artifical_vars_in_basis:
            if xB[iB.index(some_var), 0] != 0:
                return [16, iB, iN, xB]
        # Artificial variable in basis ends up with zero value, thus
        return [4, iB, iN, xB]


@protect_method
def simplex_method(A, b, c, irule):
    """
    Find a basic optimal solution for the linear program

        min:    c*x
        ST:     Ax=b
                x>=0,

    where A is an (m,n) matrix.

    Input Parameters:

    A - (n,m) constraint matrix
    b - (m,1) POSITIVE vector appearing in the constraint equation above
    c - (1,n) vector giving the coefficients of the objective function

    irule - integer parameter specifying which pivot rule to use
        - irule = 0 indicates that the smallest coefficient rule should be used
        - irule = 1 indicates that Bland's rule should be used

    Output Parameters:

    istatus - integer parameter reporting the condition of the
        - istatus = 0  indicates normal completion (i.e., a solution has been found and reported)
        - istatus = 4  indicates the program is infeasible
        - istatus = 16 indicates the program is feasible but our initialization procedure has failed
        - istatus = 32 indicates that the program is unbounded

    X  - vector of length n specifying the solution
    eta - the minimum value of the objective function
    iB - integer vector specifying the m indices of the basic variables after the simplex step
    iN - integer vector specifying the n-m indices of the nonbasic variables after the simplex step
    xB - vector of length m specifying the values of the basic variables after the simplex step
    """
    istatus = None
    X = None
    eta = 0
    iB = None
    iN = None
    xB = None

    # @@ Write program here to assign appropriate values to [istatus, X, eta, iB, iN, xB] >>>
    # Step 1: Start with first phase of simplex
    [istatus, iB, iN, xB] = simplex_init(A, b, c, irule=irule)

    # Infeasibility
    if istatus == 16:
        istatus = 4
        return [istatus, X, eta, iB, iN, xB]

    # Failure to initialize
    if istatus == 4:
        istatus = 16
        return [istatus, X, eta, iB, iN, xB]

    # Step 2: Start Phase-2
    Binv = inv(A[:, [i - 1 for i in iB]])
    phase2_termination = {-1, 16}
    while (istatus not in phase2_termination):
        [istatus, iB, iN, xB, Binv] = simplex_step(A, b, c, iB, iN, xB, Binv, irule)

    # This is unbounded case
    if istatus == 16:
        istatus = 32
        return [istatus, X, eta, iB, iN, xB]

    # This is optimal solution
    else:
        istatus = 0
        X = np.zeros((6, 1), dtype=np.float64)
        X[[(b1 - 1) for b1 in iB]] = xB
        eta = (c * X)[0, 0]

    return [istatus, X, eta, iB, iN, xB]

# Elementary Row Operations >> Rahul's code on Row Operations

# for i in range(len(A[index_dict[leaving_var_ind], :])):
#     A[index_dict[leaving_var_ind], i] /= A[index_dict[leaving_var_ind], entering_var_index - 1]
#     for j in range(len(A[:, entering_var_index - 1])):
#         for i in range(len(A[index_dict[leaving_var_ind], :])):
#             A[j, i] -= (A[index_dict[leaving_var_ind], i] * A[j, entering_var_index - 1])
#
# for i in range(c.shape[1]):
#     c[0, entering_var_index - 1] -= (
#             c[0, entering_var_index - 1] * A[index_dict[leaving_var_ind], entering_var_index - 1])
#     if i != (entering_var_index - 1):
#         c[0, i] = c[0, i] - (c[0, entering_var_index - 1] * A[index_dict[leaving_var_ind], i])
#
# b_hat[index_dict[leaving_var_ind], 0] /= A[index_dict[leaving_var_ind] - 1, entering_var_index - 1]
# for i in range(len(b)):
#     if (i != index_dict[leaving_var_ind]):
#         b_hat[i, 0] = b_hat[i, 0] - (b_hat[index_dict[leaving_var_ind], 0] * A[i, entering_var_index - 1])
# for i in range(c.shape[1]):
#     if c[0, i] < 0:
#         counter1 += 1
# if counter1 == 0:
#     istatus = -1  # optimal feasible vector found
# elif (counter1 != 0 and counter2 == 0):
#     istatus = 0  # one iteration of non-degenerate simplex

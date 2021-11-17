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

    irule - integer parameter specifying which pivot rule to use:
        irule = 0 indicates that the smallest coefficient rule should be used
        irule = 1 indicates that Bland's rule should be used

    Output Parameters:

    istatus - integer parameter reporting the condition of the
       istatus = 0  indicates normal completeion (i.e., a solution has been found and reported)
       istatus = 4  indicates the program is infeasible
       istatus = 16 indicates the program is feasible but our initialization procedure has failed
       istatus = 32 indicates that the program is unbounded

    X  - vector of length n specifying the solution
    eta - the minimum value of the objective function
    iB - integer vector specifying the m indices of the basic variables after the simplex step
    iN - integer vector specifying the n-m indices of the nonbasic variables after the simplex step
    xB - vector of length m specifying the values of the basic variables after the simplex step
    """
    return


def simplex_init(A, b, c):
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
    pass


def simplex_step(A, b, c, iB, iN, xB, irule):
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
    return [A, b, c, iB, iN, xB, irule]

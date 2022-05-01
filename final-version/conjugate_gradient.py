import numpy as np

def conjugate_gradient(A, b, n):
    '''
    Serial conjugate gradient solver (Dense Matrices).
    Input: 
        A - nxn numpy array (variable coefficients)
        b - nx1 numpy array (y-intercept)
        n - number of iterations
    Return:
        x - nx1 numpy array solution
    '''

    N = n**2
    x = np.zeros(N)
    r = b - np.dot(A, x)
    p = r
    rsold = np.dot(r, r)
    for _ in range(N):
        z = np.dot(A, p)
        alpha = rsold / np.dot(p, z)
        x = x + alpha * p
        r = r - alpha * z
        rsnew = np.dot(r, r)
        if rsnew < 1e-10:
            break
        p = r + (rsnew / rsold) * p
        rsold = rsnew
    
    return x

def sanity_check(A, b, n, x, sol):
    '''
    Sanity check for the conjugate gradient solver.
    '''
    assert np.all(sol == x)


    print(
'''
-----------------------------------------------------
Inputs:
1. variable coefficients - 
{}
2. y-intercept - 
{}
3.number of iterations - 
{}
Output:
Provided Output - 
{}

Validity:
1. Expected Output - 
{}
2. Program Success - 
{}
-----------------------------------------------------
'''
.format(A, b.T, n, x.T, sol.T, np.all(sol == x)))

    

if __name__ == '__main__':
    n = 3
    N = n**2
    A = np.array([np.ones(N) for _ in np.arange(N)]).T
    b = np.ones(N)

    x = conjugate_gradient(A, b, n)
    sol = np.array([1/N for _ in np.arange(N)])
    
    sanity_check(A, b, n, x, sol)


    


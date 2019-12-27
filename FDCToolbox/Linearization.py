from F16_Model import f16_model
from Trim import trim
from utils import print_matrix
import numpy as np


# Ref [1] page 203-204
# Ref [3] https://github.com/stanleybak/AeroBenchVVPython/blob/master/code/jacobFun.py
def numjac(Xequil, Uequil, Xcg, printOn=False):
    '''
    Numerical Jacobian Function
    numerically calculates the linearized A, B, C, & D matrices
    of an F-16 model given certain parameters where:

    x = A x  +  B u
    y = C x  +  D u

    There are 13 state variables, 4 inputs, and 10 outputs, making the matrix sizes:
    A: 13x13, B: 13x4, C: 10x13, D: 10x4
    '''

    if printOn:
        print('Running jacobFun.py')

    xe = Xequil.copy()
    ue = Uequil.copy()
    x = xe.copy()
    u = ue.copy()

    n = len(Xequil)
    m = len(Uequil) 

    tol = 1e-6

    xde = f16_model(x, u, Xcg=Xcg)

    #####   A matrix    #####
    dx = 0.01*x
    for i in range(0,n):
        if dx[i] == 0.0:
            dx[i] = 0.1

    last = np.zeros((n,1), dtype=float)
    A = np.zeros((n,n), dtype=float)

    for j in range(0,n):
        xt = x
        for i in range(0,10):
            xt[j] = x[j] + dx[j]
            xd1 = f16_model(xt, u, Xcg=Xcg)
            xt[j] = x[j] - dx[j]
            xd2 = f16_model(xt, u, Xcg=Xcg)
            A[:, j] = (np.transpose(xd1.ravel() - xd2.ravel()) / (2*dx[j]))
            if np.max(np.abs(A[:,j] - last) / abs(A[:,j] + 1e-12)) < tol:
                break
            dx[j] = 0.5*dx[j]
            last = A[:,j]
        ## column = j
        iteration = i
        if iteration == 10:
            print(f"not converged on A, column {j}")

    AA = 2*A

    #####   B matrix    #####
    du = 0.01*u

    for i in range(0,m):
        if du[i] == 0.0:
            du[i] = 0.1

    last = np.zeros((n,1), dtype=float)
    B = np.zeros((n,m), dtype=float)

    for j in range(0,m):
        ut = u
        for i in range(0,10):
            ut[j] = u[j] + du[j]
            ut = ut
            xd1 = f16_model(x, ut, Xcg=Xcg)
            ut[j] = u[j] - du[j]
            ut = ut
            xd2 = f16_model(x, ut, Xcg=Xcg)
            B[:, j] = (np.transpose(xd1.ravel() - xd2.ravel()) / (2*du[j]))
            if np.max(np.abs(B[:,j] - last) / abs(B[:,j] + 1e-12)) < tol:
                break
            dx[j] = 0.5*dx[j]
            du[j] = 0.5*du[j]
            last = B[:,j]
        ## column = j
        iteration = i
        if iteration == 10:
            print(f"not converged on B, column {j}")    
    
    BB = 2*B

    return AA, BB


def linearize(orient, inputs, Xcg, printOn=False):
    """
    Get linearized version of f16 model about a trim point
    """
    
    Xequil, Uequil = trim(orient, inputs, Xcg, printOn)

    A, B = numjac(Xequil, Uequil, Xcg, printOn)        

    if printOn:
        print("A_matrix = ")
        print_matrix(A)
        print(" ")
        print("B_matrix = ")
        print_matrix(B)
        print(" ")
    
    return A, B


def get_lon_A(A_matrix):
    
    A = A_matrix

    A_lon = np.array(
        [
        [A[0,0], A[0,1], A[0,4], A[0,7]],
        [A[1,0], A[1,1], A[1,4], A[1,7]],
        [A[4,0], A[4,1], A[4,4], A[4,7]],
        [A[7,0], A[7,1], A[7,4], A[7,7]]
        ]
    )

    return A_lon

def get_lon_B(B_matrix):

    B = B_matrix

    B_lon = np.array(
        [
        [B[0,0], B[0,1]],
        [B[1,0], B[1,1]],
        [B[4,0], B[4,1]],
        [B[7,0], B[7,1]]
        ]
    )

    return B_lon

def get_lat_A(A_matrix):

    A = A_matrix

    A_lat = np.array(
        [
        [A[2,2], A[2,3], A[2,6], A[2,8]],
        [A[3,2], A[3,3], A[3,6], A[3,8]],
        [A[6,2], A[6,3], A[6,6], A[6,8]],
        [A[8,2], A[8,3], A[8,6], A[8,8]]
        ]
    )

    return A_lat

def get_lat_B(B_matrix):

    B = B_matrix

    B_lat = np.array(
        [
        [B[2,2], B[2,3]],
        [B[3,2], B[3,3]],
        [B[6,2], B[6,3]],
        [B[8,2], B[8,3]]
        ]
    )

    return B_lat


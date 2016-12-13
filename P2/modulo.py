from numpy import *

"""
Computes the modulo of all elements in the matrix A with the positive number n. If an element is negative then the
modulo will also be negative.
"""
def mod(A,n):
    B = A
    for i in range(0,shape(A)[0]):
        for j in range(0,shape(A)[1]):
            if A[i,j] < 0:
                B[i,j] = -(-A[i,j]%n)
            else:
                B[i,j] = A[i,j]%n
    return B


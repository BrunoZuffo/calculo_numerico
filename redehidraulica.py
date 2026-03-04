import numpy as np

conec = np.array([
    [1,2],
    [2,3],
    [3,4],
    [4,5],
    [5,2],
    [5,3],
    [5,1]
])

C = np.array([2.0, 2.0, 1.0, 2.0, 1.0, 2.0, 2.0])

def Assembly(conec, C):
    nv = 5
    nc = 7
    A = np.zeros(shape=(nv,nv))
    for k in range(nc):
        n1 = conec[k,0] - 1
        n2 = conec[k,1] - 1

        A[n1,n1] += C[k]
        A[n2,n2] += C[k]
        A[n1,n2] -= C[k]
        A[n2,n1] -= C[k]

    return A

print(Assembly(conec, C))
    
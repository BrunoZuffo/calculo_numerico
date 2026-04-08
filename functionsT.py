import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union
from scipy import sparse
from scipy.sparse.linalg import spsolve

def ij2n (i, j, N):
    return i + j*N

# FUNÇÃO ASSEMBLY --------------------------------------------------------------------------
def Assembly(N, k):
    nunk = N * N # Número total de pontos/incógnitas 
    A = np.zeros(shape=(nunk, nunk)) # Cria a matriz inicial cheia de zeros 
    
    # Varre apenas os nós internos da malha 
    for i in range(1, N-1):
        for j in range(1, N-1):
            
            # 1. Calcula o índice global do ponto central (Ic) 
            Ic = ij2n(i, j, N)
            
            # 2. Calcula os índices globais dos 4 vizinhos 
            Ie = ij2n(i+1, j, N) # Vizinho Leste (East)
            Iw = ij2n(i-1, j, N) # Vizinho Oeste (West)
            In = ij2n(i, j+1, N) # Vizinho Norte (North)
            Is = ij2n(i, j-1, N) # Vizinho Sul (South)
            
            # 3. Preenche a matriz A na linha correspondente ao ponto central (Ic) 
            # O ponto central recebe 4*k, e os vizinhos recebem -k 
            A[Ic, [Ic, Ie, Iw, In, Is]] = [4*k, -k, -k, -k, -k]
            
    return A

# FUNÇÃO SolveSystem --------------------------------------------------------------------------
def SolveSystem(N, h, k, TL, TR, TB, TT, fonte):
    
    A = Assembly(N, k)
    Atilde = A.copy()
    
    nunk = N**2
    b = np.zeros(nunk)
    
    if fonte is not None:
        b += fonte * h**2
    
    for i in range(N):
        for j in range(N):
            Ic = ij2n(i, j, N)
            
            if i == 0:  # esquerda
                Atilde[Ic, :] = 0
                Atilde[Ic, Ic] = 1
                b[Ic] = TL
                
            elif i == N-1:  # direita
                Atilde[Ic, :] = 0
                Atilde[Ic, Ic] = 1
                b[Ic] = TR
                
            elif j == 0:  # base
                Atilde[Ic, :] = 0
                Atilde[Ic, Ic] = 1
                b[Ic] = TB[i]   # usa vetor
                
            elif j == N-1:  # topo
                Atilde[Ic, :] = 0
                Atilde[Ic, Ic] = 1
                b[Ic] = TT[i]   # usa vetor
    
    T = np.linalg.solve(Atilde, b)
    
    return T.reshape((N, N))

# FUNÇÃO SolveSystemSparse --------------------------------------------------------------------------
def SolveSystemSparse(N, h, k, TL, TR, TB, TT, fonte):
    
    nunk = N * N
    
    # 1. Construção da matriz A
    d1 = np.ones(nunk) * 4.0 * k
    d2 = -np.ones(nunk - 1) * k
    d3 = -np.ones(nunk - N) * k

    A = sparse.diags(
        [d3, d2, d1, d2, d3],
        [-N, -1, 0, 1, N],
        format='lil'
    )
    
    # 2. Vetor b
    b = np.zeros(nunk)
    
    if fonte is not None:
        b += fonte * h**2
    
    # 3. Condições de contorno (ATUALIZADO)
    for i in range(N):
        for j in range(N):
            Ic = ij2n(i, j, N)
            
            if i == 0:  # esquerda
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TL
                
            elif i == N-1:  # direita
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TR
                
            elif j == 0:  # base
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TB[i]   # vetor
                
            elif j == N-1:  # topo
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TT[i]   # vetor
    
    # 4. Converter para CSR
    A = A.tocsr()
    
    # 5. Resolver
    T = spsolve(A, b)
    
    return T.reshape((N, N))

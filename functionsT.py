import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union

def ij2n (i, j, N):
    return i + j*N


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
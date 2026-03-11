import matplotlib
matplotlib.use("TkAgg")

from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, createK, createD, calc_vazao, calc_potencia
import numpy as np
import matplotlib.pyplot as plt


# matriz de conectividade do exemplo
conec = np.array([
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 2],
    [5, 3],
    [5, 1]
])


C = np.array([2, 2, 1, 2, 1, 2, 2]) # valores de Ck do exemplo

Xno = np.array([
    [0, 0],   # nó 1
    [1, 0],   # nó 2
    [2, 0],   # nó 3
    [3, 0],   # nó 4
    [1.5, 1]  # nó 5
])

matriz = Assembly(conec, C)

print(matriz)

pressure = SolveNetwork(conec,C,3,1,3)

print(pressure)

matriz_vazao = calc_vazao(conec, C, pressure)
print (matriz_vazao)
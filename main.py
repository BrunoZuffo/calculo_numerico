from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, calc_vazao
import numpy as np

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

matriz = Assembly(conec, C)

print(matriz)

pressure = SolveNetwork(conec,C,3,1,3)

print(pressure)

matriz_vazao = calc_vazao(conec, C, pressure)
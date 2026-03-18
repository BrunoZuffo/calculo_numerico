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

#teste para exercicio 2
ps = {
    1: 100,   # nó 1 → alta pressão (entrada)
    4: 0      # nó 4 → baixa pressão (saída)
}

Qs = {
    3: 10     # nó 3 → injeção de vazão
}

matriz = Assembly(conec, C)

print("MATRIZ A")
print(matriz)

pressure = SolveNetwork(conec,C, ps=ps, Qs=Qs)

print("PRESSURE:")
print(pressure)

matriz_vazao = calc_vazao(conec, C, pressure)
print("MATRIZ DE VAZÃO:")
print (matriz_vazao)

fig, ax = PlotaRede(conec, Xno, pressure, matriz_vazao)
plt.show
from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork
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

matriz = SolveNetwork(conec,C,3,1,3)

print(matriz)

def SolveNetworkAtualizada(Conec,C,natm,Qs):
    natm = natm - 1

    Atilde = Assembly(conec, C)

    # condição de pressão atmosférica
    Atilde[natm, :] = 0
    Atilde[natm, natm] = 1

    n = Atilde.shape[0]

    b = np.zeros(n)

    # adicionar todas as vazões injetadas
    for node, q in Qs.items():
        b[node] = q

    pressure = np.linalg.solve(Atilde, b)

    return pressure
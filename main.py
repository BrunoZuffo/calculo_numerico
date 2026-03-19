import matplotlib
matplotlib.use("TkAgg")

from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, createK, createD, calc_vazao, calc_potencia, AssemblyVectorC
import numpy as np
import matplotlib.pyplot as plt

## matriz de conectividade do exemplo
#conec = np.array([
#    [1, 2],
#    [2, 3],
#    [3, 4],
#    [4, 5],
#    [5, 2],
#    [5, 3],
#    [5, 1]
#])
#
##coordenadas dos nós
#Xno = np.array([
#    [0, 0],   # nó 1
#    [1, 0],   # nó 2
#    [2, 0],   # nó 3
#    [3, 0],   # nó 4
#    [1.5, 1]  # nó 5
#])
#
## valores de Ck do exemplo
#C = np.array([2, 2, 1, 2, 1, 2, 2]) 

#teste para exercicio 2
ps = {
    1: 100,   # nó 1 → alta pressão (entrada)
    4: 0      # nó 4 → baixa pressão (saída)
}

Qs = {
    3: 10,     # nó 3 → injeção de vazão
    100: 200
}

Xno, conec = GeraGrafo(levels=3)

mm_to_m = 0.001
Xno = Xno * mm_to_m

C = AssemblyVectorC(conec, Xno)

natm = len(Xno) - 1  #nó atmosférico
nbomba = 0  #nó conectado à bomba
Qbomba = 1.0e-7  #vazao da bomba

Qs = {nbomba: Qbomba} # dicionario

matriz = Assembly(conec, C)

print(matriz)

pressure = SolveNetwork(conec,C,natm,Qs)

print(pressure)

matriz_vazao = calc_vazao(conec, C, pressure)
print (matriz_vazao)

fig, ax = PlotaRede(conec, Xno, pressure, matriz_vazao, factor_units=mm_to_m)
plt.show()

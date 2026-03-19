import matplotlib
matplotlib.use("TkAgg")

from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, createK, createD, calc_vazao, calc_potencia, AssemblyVectorC
import numpy as np
import matplotlib.pyplot as plt

Xno, conec = GeraGrafo(levels=10)

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

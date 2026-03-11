from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork
import numpy as np

C = np.array([2, 2, 1, 2, 1, 2, 2]) # valores de Ck do exemplo

<<<<<<< HEAD

#gerar grafo
Xno, conec = GeraGrafo(levels=3)
=======
matriz = Assembly(conec, C)
>>>>>>> main

#converter as unidades
mm_to_m = 0.001
Xno = Xno * mm_to_m

<<<<<<< HEAD
#criando condutancia para cada tubo
C = np.ones(len(conec))

#resolver sistema - gera grafo devolve conec na base 0, mas dentro da assembly e solvenetwork o temos conec-1
pressure = SolveNetwork(conec + 1, C, 3, 1, 3)

#calcular vazao
q = calc_vazao(conec + 1, C, pressure)

fig, ax = PlotaRede(conec, Xno, pressure, q, mm_to_m)
plt.show()
=======
matriz = SolveNetwork(conec,C,3,1,3)

print(matriz)
>>>>>>> main

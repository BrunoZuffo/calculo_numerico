from functionsT import ij2n, Assembly, SolveSystem, SolveSystemSparse
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import time

# Valores de teste:

Lx = 0.02  # 2 cm
Ly = 0.01  # 1 cm

N = 50
h = Lx / (N - 1)

K = 1

TL = 10
TR = 30

x_coords = np.linspace(0, Lx, N)

TB = 10 + 20 * (x_coords / Lx)
TT = 10 + 20 * (x_coords / Lx)

fonte = None

# Teste Função Assembly
A = Assembly(N, K)
print("FUNÇÃO ASSEMBLY:")
print(A)
print("")

# Teste Função SolveSystem
print("FUNÇÃO SolveSystem:")
T = SolveSystem(N, h, K, TL, TR, TB, TT, fonte)
print(T)
print("")

# Teste Função SolveSystemSparse
print("FUNÇÃO SolveSystemSparse:")
T1 = SolveSystemSparse(N, h, K, TL, TR, TB, TT, fonte)
print(T1)
print("")
from functionsT import ij2n, Assembly, SolveSystem, SolveSystemSparse, PlotaPlaca, Jacobi, GaussSeidel
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import time

# INÍCIO ----------------------------------------------------------------------------------
# não faz parte de nenhum exercício

# Valores da placa para casos de teste:

Nx = 51
Ny = 26

# Valores da placa assumidos, considerando o material sendo o polidimetilsiloxano puro (PDMS):
# esses valores são fixos

Lx = 0.02  # 2 cm
Ly = 0.01  # 1 cm

h = Lx / (Nx - 1)

K = 0.25 #valor intermediário do intervalo: 0.2 − 0.3 W ·K^−1 · m^−1

TL = 10 # graus Celsius
TR = 30 # graus Celsius

fonte = 5.0e5  # valor do intervalo: 10^5 – 10^6 W ·m^−3 - NÃO ALTERAR

x_coords = np.linspace(0, Lx, Nx)

TB = 10 + 20 * (x_coords / Lx)
TT = 10 + 20 * (x_coords / Lx)

# Teste Função Assembly
A = Assembly(Nx, Ny, K)
print("FUNÇÃO ASSEMBLY:")
print(A)
print("")

# Teste Função SolveSystem
#print("FUNÇÃO SolveSystem:")
#T1, t_assembly, t_montagem, t_sistema = SolveSystem(Nx, Ny, h, K, TL, TR, TB, TT, fonte)
#print(T)
#print("")

# Teste Função SolveSystemSparse
print("FUNÇÃO SolveSystemSparse:")
T, t_assembly, t_montagem, t_sistema= SolveSystemSparse(Nx, Ny, h, K, TL, TR, TB, TT, fonte)
print(T)
print("")

# Teste da Função PlotaPlaca
PlotaPlaca(Nx=Nx, Ny=Ny, Lx=Lx, Ly=Ly, T=T, flag_type='contour')

# Exercício 1 -------------------------------------------------------------------------------

casos = [(21,11), (41,21), (81,41), (161,81), (321,161)]
resultados = []

for (Nx, Ny) in casos:

    h = Lx / (Nx - 1)
    x_coords = np.linspace(0, Lx, Nx)
    TB = 10 + 20 * (x_coords / Lx)
    TT = 10 + 20 * (x_coords / Lx)
    
    # DENSO
    if Nx <= 161: # o computador não tem memória suficiente para resolver o caso 321 X 161 sem usar matriz esparsa
        T_d, tA_d, tM_d, tS_d = SolveSystem(Nx, Ny, h, K, TL, TR, TB, TT, fonte)
    else:
        T_d, tA_d, tM_d, tS_d = None, None, None, None
    
    # ESPARSO
    T_s, tA_s, tM_s, tS_s = SolveSystemSparse(Nx, Ny, h, K, TL, TR, TB, TT, fonte)
    
    resultados.append([
        Nx, Ny,
        tA_d, tM_d, tS_d,
        tA_s, tM_s, tS_s
    ])
    
    # CONTOUR
    PlotaPlaca(Nx, Ny, Lx, Ly, T_s)
    
    # PERFIL CENTRAL (gráfio da temperatura em função do eixo X)
    linha_central = Ny // 2
    perfil = T_s[linha_central, :]
    
    x = np.linspace(0, Lx, Nx)
    
    plt.plot(x, perfil)
    plt.title(f'Temperatura ao Longo do Eixo Central ({Nx} X {Ny})')
    plt.xlabel('x (m)')
    plt.ylabel('Temperatura (°C)')
    plt.grid()
    plt.show()

print("\n" + "="*90)
print(f"{'COMPARAÇÃO DE TEMPOS (DENSO vs ESPARSO)':^90}")
print("="*90)

header = (
    f"{'Nx':>5} {'Ny':>5} | "
    f"{'Mont.(D)':>10} {'Sist.(D)':>10} {'Resol.(D)':>10} | "
    f"{'Mont.(S)':>10} {'Sist.(S)':>10} {'Resol.(S)':>10}"
)

print(header)
print("-"*90)

for r in resultados:
    print(
        f"{r[0]:5d} {r[1]:5d} | "
        f"{(f'{r[2]:.4f}' if r[2] is not None else '---'):>10} "
        f"{(f'{r[3]:.4f}' if r[3] is not None else '---'):>10} "
        f"{(f'{r[4]:.4f}' if r[4] is not None else '---'):>10} | "
        f"{(f'{r[5]:.4f}' if r[5] is not None else '---'):>10} "
        f"{(f'{r[6]:.4f}' if r[6] is not None else '---'):>10} "
        f"{(f'{r[7]:.4f}' if r[7] is not None else '---'):>10}"
    )

print("="*90)

# Exercício 1, parte 2 -------------------------------------------------------------------------------

casos = [(21,11), (41,21), (81,41), (161,81), (321,161)]
TOLs = [1e-2, 1e-4, 1e-6]
MAXITER = 10000

print("\n" + "="*90)
print(f"{'COMPARAÇÃO DE MÉTODOS (JACOBI vs GAUSS-SEIDEL)':^90}")
print("="*90)

print(f"{'Nx':>5} {'Ny':>5} {'TOL':>8} | {'Jacobi(s)':>10} {'Iter':>6} | {'GS(s)':>10} {'Iter':>6}")

for (Nx, Ny) in casos:
    
    x_coords = np.linspace(0, Lx, Nx)
    TB = 10 + 20 * (x_coords / Lx)
    TT = 10 + 20 * (x_coords / Lx)
    h = Lx / (Nx - 1)
    
    for tol in TOLs:
        
        Tj, it_j, tj = Jacobi(Nx, Ny, h, K, TL, TR, TB, TT, fonte, tol, MAXITER)
        Tg, it_g, tg = GaussSeidel(Nx, Ny, h, K, TL, TR, TB, TT, fonte, tol, MAXITER)
        
        print(f"{Nx:5d} {Ny:5d} {tol:8.0e} | {tj:10.4f} {it_j:6d} | {tg:10.4f} {it_g:6d}")
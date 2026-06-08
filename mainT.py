

from functionsT import ij2n, Assembly, SolveSystem, SolveSystemSparse, PlotaPlaca , SolveSystemSparse_Circle
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import time

# INÍCIO ----------------------------------------------------------------------------------
# não faz parte de nenhum exercício

# Valores da placa para casos de teste:

Nx = 50
Ny = 20

# Valores da placa assumidos, considerando o material sendo o polidimetilsiloxano puro (PDMS):
# esses valores são fixos

Lx = 0.02  # 2 cm
Ly = 0.01  # 1 cm

h = Lx / (Nx - 1)

K = 0.25 #valor intermediário do intervalo: 0.2 − 0.3 W ·K^−1 · m^−1

TL = 10 # graus Celsius
TR = 30 # graus Celsius

fonte = 5.0e5  # valor do intervalo: 10^5 – 10^6 W ·m^−3

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


# Exercício 5 ----------------------------------------------------------------------------------------------

print("\n" + "=" * 75)
print(f"{'RESULTADOS DO EX5  (Coeficientes a, b e c)':^75}")
print("=" * 75)

# 1. Definir a malha para o Exercício 5 (conforme sugerido no material)
Nx_5, Ny_5 = 101, 51 
h_5 = Lx / (Nx_5 - 1)
x_coords_5 = np.linspace(0, Lx, Nx_5)
TB_5 = 10 + 20 * (x_coords_5 / Lx)
TT_5 = 10 + 20 * (x_coords_5 / Lx)

# --- CORREÇÃO: Parâmetros da região circular declarados aqui ---
R  = 0.002          # raio de 0.2 cm
xc = 0.75 * Lx      # posição X do centro
yc = 0.50 * Ly      # posição Y do centro
# ---------------------------------------------------------------

# 2. Escolher um ponto k interno arbitrário (longe do círculo e das bordas)
# Vamos pegar o meio da altura (Ny//2) e um quarto do comprimento (Nx//4)
i_k = Nx_5 // 4
j_k = Ny_5 // 2

# CENÁRIO 1: TR e TC com valores nominais (30°C)
TR_1, TC_1 = 30.0, 30.0
T_s1, _, _, _, _ = SolveSystemSparse_Circle(
    Nx_5, Ny_5, h_5, K, TL, TR_1, TB_5, TT_5, fonte,
    Lx, Ly, R, xc, yc, TC_1
)
Tk_1 = T_s1[j_k, i_k]

# CENÁRIO 2: Aumentar apenas TR (ex: para 40°C)
TR_2, TC_2 = 40.0, 30.0
T_s2, _, _, _, _ = SolveSystemSparse_Circle(
    Nx_5, Ny_5, h_5, K, TL, TR_2, TB_5, TT_5, fonte,
    Lx, Ly, R, xc, yc, TC_2
)
Tk_2 = T_s2[j_k, i_k]

# CENÁRIO 3: Aumentar apenas TC (ex: para 40°C)
TR_3, TC_3 = 30.0, 40.0
T_s3, _, _, _, _ = SolveSystemSparse_Circle(
    Nx_5, Ny_5, h_5, K, TL, TR_3, TB_5, TT_5, fonte,
    Lx, Ly, R, xc, yc, TC_3
)
Tk_3 = T_s3[j_k, i_k]

# CÁLCULO DOS COEFICIENTES (Equação: Tk = a*TR + b*TC + c)
# a = Variação de Tk / Variação de TR
a = (Tk_2 - Tk_1) / (TR_2 - TR_1)

# b = Variação de Tk / Variação de TC
b = (Tk_3 - Tk_1) / (TC_3 - TC_1)

# c = Isolando a constante no Cenário 1
c = Tk_1 - (a * TR_1) - (b * TC_1)

print(f"Ponto k escolhido: índice [i={i_k}, j={j_k}] (x = {i_k*h_5*100:.2f} cm, y = {j_k*h_5*100:.2f} cm)")
print(f"Temperatura Tk no Cenário Nominal (TR=30, TC=30): {Tk_1:.4f} °C\n")
print("Coeficientes encontrados da relação [Tk = a*TR + b*TC + c]:")
print(f"a = {a:.6f}")
print(f"b = {b:.6f}")
print(f"c = {c:.6f}\n")
print(f"Equação Final: Tk = {a:.4f}*TR + {b:.4f}*TC + {c:.4f}")
print("=" * 75)
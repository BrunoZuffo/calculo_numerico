from functionsT import ij2n, Assembly, SolveSystem, SolveSystemSparse, SolveSystemSparse_Circle, PlotaPlaca, Jacobi, GaussSeidel, AnimacaoTemperatura
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
#print("FUNÇÃO ASSEMBLY:")
#print(A)
#print("")

# Teste Função SolveSystem
#print("FUNÇÃO SolveSystem:")
#T1, t_assembly, t_montagem, t_sistema = SolveSystem(Nx, Ny, h, K, TL, TR, TB, TT, fonte)
#print(T)
#print("")

# Teste Função SolveSystemSparse
#print("FUNÇÃO SolveSystemSparse:")
T, t_assembly, t_montagem, t_sistema= SolveSystemSparse(Nx, Ny, h, K, TL, TR, TB, TT, fonte)
#print(T)
#print("")

# Teste da Função PlotaPlaca
PlotaPlaca(Nx=Nx, Ny=Ny, Lx=Lx, Ly=Ly, T=T, flag_type='contour')

# 2.5.1 Exercício 1 -------------------------------------------------------------------------------

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

# 2.5.1 Exercício 2 ----------------------------------------------------------------------------------------------

# Parâmetros da região circular
R  = 0.002          
xc = 0.75 * Lx     
yc = 0.50 * Ly     
TC = 30.0           

casos_ex2      = [(21, 11), (41, 21), (81, 41), (161, 81), (321, 161)]
resultados_ex2 = []
 
for (Nx, Ny) in casos_ex2:
 
    h        = Lx / (Nx - 1)
    x_coords = np.linspace(0, Lx, Nx)
    TB       = 10 + 20 * (x_coords / Lx)
    TT       = 10 + 20 * (x_coords / Lx)
 
    T_s, tA_s, tM_s, tS_s, circle_mask = SolveSystemSparse_Circle(
        Nx, Ny, h, K, TL, TR, TB, TT, fonte,
        Lx, Ly, R, xc, yc, TC
    )
 
    T_max = T_s.max()
    resultados_ex2.append([Nx, Ny, tA_s, tM_s, tS_s, T_max])
 
    # Curvas de nível 
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
 
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_aspect('equal')
    levels = np.linspace(10, 42, 17)
    im  = ax.contourf(X, Y, T_s, levels=levels, cmap='jet')
    ax.contour(X, Y, T_s, levels=levels, linewidths=0.25, colors='k')
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.22)
    cbar.set_ticks(np.arange(10, 46, 4))
    ax.set(xlabel='x', ylabel='y', title=f'Contours of temperature  ({Nx}×{Ny})')
    ax.set_xticks([0, 0.01, 0.02])
    ax.set_yticks([0.000, 0.005, 0.010])
    plt.tight_layout()
    plt.show()

    
    # Perfil de temperatura ao longo do eixo central
    linha_central = Ny // 2
    perfil        = T_s[linha_central, :]
 
    plt.figure(figsize=(7, 3.5))
    plt.plot(x * 100, perfil, color='crimson', lw=1.8)
    plt.axvline(xc * 100, color='steelblue', ls='--', lw=1.2,
                label=f'Centro do círculo  (x = {xc*100:.1f} cm)')
    plt.title(f'Perfil de temperatura – eixo central  ({Nx}×{Ny})')
    plt.xlabel('x (cm)')
    plt.ylabel('Temperatura (°C)')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()
 
# Tabela de resultados 
print("\n" + "=" * 90)
print(f"{'RESULTADOS DO EX2  (Região Circular  TC = 30°C)':^90}")
print("=" * 90)
 
header2 = (
    f"{'Nx':>5} {'Ny':>5} | "
    f"{'Mont. (S)':>10} {'Sist.(S)':>10} {'Resol. (sS)':>10} | "
    f"{'T_max (°C)':>10}"
)
print(header2)
print("-" * 90)
 
for r in resultados_ex2:
    print(
        f"{r[0]:5d} {r[1]:5d} | "
        f"{r[2]:10.4f} {r[3]:10.4f} {r[4]:10.4f} | "
        f"{r[5]:10.3f}"
    )
 
print("=" * 90)

# 2.5.1 Exercício 4 ----------------------------------------------------------------------------------------------

# Definindo os parâmetros da malha fixa 
Nx_param, Ny_param = 101, 51
h_param = Lx / (Nx_param - 1)

x_coords_param = np.linspace(0, Lx, Nx_param)
TB_param = 10 + 20 * (x_coords_param / Lx)
TT_param = 10 + 20 * (x_coords_param / Lx)

# Vetores de Temperatura do Círculo (TC) para testar
valores_TC = np.linspace(0, 100, 21) 

T_max_list = []
T_med_list = []

# Loop para rodar a simulação para cada valor de TC
for TC_atual in valores_TC:
    
    T_vec, _, _, _, _ = SolveSystemSparse_Circle(
        Nx_param, Ny_param, h_param, K, TL, TR, TB_param, TT_param, fonte, Lx, Ly, R, xc, yc, TC_atual
    )
    
    # O SciPy retorna o vetor 1D. Achamos o máximo e a média com o NumPy
    T_max_list.append(np.max(T_vec))
    T_med_list.append(np.mean(T_vec))

# Geração dos Gráficos
plt.figure(figsize=(8, 5))

# Eixo Y esquerdo (Temperatura Máxima)
ax1 = plt.gca()
linha1, = ax1.plot(valores_TC, T_max_list, color='crimson', marker='o', label='Temp. Máxima')
ax1.set_xlabel('Temperatura do Círculo TC (°C)')
ax1.set_ylabel('Temperatura Máxima da Placa (°C)', color='crimson')
ax1.tick_params(axis='y', labelcolor='crimson')

# Eixo Y direito (Temperatura Média)
ax2 = ax1.twinx() 
linha2, = ax2.plot(valores_TC, T_med_list, color='steelblue', marker='s', label='Temp. Média')
ax2.set_ylabel('Temperatura Média da Placa (°C)', color='steelblue')
ax2.tick_params(axis='y', labelcolor='steelblue')

# Título e Legendas
plt.title('Estudo Paramétrico: Efeito de TC na Placa (Malha 101x51)')
ax1.grid(True, linestyle='--', alpha=0.6)

# Juntar as legendas dos dois eixos
linhas = [linha1, linha2]
labels = [l.get_label() for l in linhas]
ax1.legend(linhas, labels, loc='upper left')

plt.tight_layout()
plt.show()

# 2.5.1 Exercício 5 ----------------------------------------------------------------------------------------------

print("\n" + "=" * 90)
print(f"{'RESULTADOS DO EX5  (Coeficientes a, b e c)':^90}")
print("=" * 90)

# 1. Definir a malha para o Exercício 5 (conforme sugerido no material)
Nx_5, Ny_5 = 101, 51 
h_5 = Lx / (Nx_5 - 1)
x_coords_5 = np.linspace(0, Lx, Nx_5)
TB_5 = 10 + 20 * (x_coords_5 / Lx)
TT_5 = 10 + 20 * (x_coords_5 / Lx)

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
print("=" * 90)

# 2.6.1 Exercício 1 -------------------------------------------------------------------------------

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
        
        Tj, it_j, tj, frame = Jacobi(Nx, Ny, h, K, TL, TR, TB, TT, fonte, tol, MAXITER, animation = False, frame_skip = 0)
        Tg, it_g, tg, frame = GaussSeidel(Nx, Ny, h, K, TL, TR, TB, TT, fonte, tol, MAXITER, animation = False, frame_skip = 0)
        
        print(f"{Nx:5d} {Ny:5d} {tol:8.0e} | {tj:10.4f} {it_j:6d} | {tg:10.4f} {it_g:6d}")

# 2.6.1 Exercício 2 -------------------------------------------------------------------------------

Nx = 41 # quando o matheus mexer nas funções de Jacobi e Gaus
Ny = 21 # mudar para 101 X 51

Lx = 0.02  # 2 cm
Ly = 0.01  # 1 cm

h = Lx / (Nx - 1)
k = 0.25

TL = 10
TR = 30

# contornos variáveis
x_coords = np.linspace(0, Lx, Nx)
TB = 10 + 20 * (x_coords / Lx)
TT = 10 + 20 * (x_coords / Lx)

fonte = 5.0e5

T, it, _, frames = GaussSeidel(Nx, Ny, h, k, TL, TR, TB, TT, fonte, TOL=1e-6, MAXIT=10000, animation=True, frame_skip = 20)
AnimacaoTemperatura(frames, Nx, Ny, Lx, Ly)

T, it, _, frames = Jacobi(Nx, Ny, h, k, TL, TR, TB, TT, fonte, TOL=1e-6, MAXIT=10000, animation=True, frame_skip = 20)
AnimacaoTemperatura(frames, Nx, Ny, Lx, Ly)
import numpy as np
import matplotlib.pyplot as plt
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.interpolate import RegularGridInterpolator
import time

from functionsHT import (ObterCondutividadeFaces_ViaNos, CriarSistemaSolidoCondutividadeVariavel,
                         calcular_temperatura_media_arestas, atualiza_condutancias, 
                         plot_nos_cromaticos_hidraulica, plot_arestas_cromaticas_hidraulics, 
                         calcular_viscosidade)

from functions import GeraGrafo, Assembly as AssemblyHidraulico

# -----------------------------------------------------------------------------
# CÓDIGO BASE - define os parâmetros da placa e executa as funções base
# -----------------------------------------------------------------------------

# Parâmetros geométricos da placa
Lx, Ly = 0.03, 0.015
k_0 = 0.25
fonte_calor = 5.0e5
TL, TR = 10.0, 30.0

# Parâmetros da inclusão circular
R_incl = 0.0025
xincl, yincl = 0.02 + R_incl, 0.0075
TC = 35.0

# Parâmetros geométricos e operacionais dos microcanais
Area_canal = 2.5e-7  # 500 um x 500 um
Q_in_0 = 1.0e-7      # 0.1 mL/s (Inlet nó 0)
Q_in_175 = 1.0e-6    # 1.0 mL/s (Inlet nó 175)

# Construção da topologia da rede hidráulica fractal
Xno, conec = GeraGrafo(levels=3)
Xno = Xno * 0.001       # Conversão de mm para metros
Xno[:, 1] += 0.5 * Ly   # Centralização geométrica no eixo Y da placa

# Nós de controle de contorno hidráulico
n_inlet_0 = 0
n_inlet_175 = 175
n_outlet = 5 

# Resolução da hidráulica nominal invariante (Caso Isotérmico de Referência a T = 20°C)
mu_iso = calcular_viscosidade(20.0)
D_h = np.sqrt(4 * Area_canal / np.pi)

C_iso = np.zeros(len(conec))
for k in range(len(conec)):
    Lk = np.linalg.norm(Xno[int(conec[k,1])] - Xno[int(conec[k,0])])
    C_iso[k] = (np.pi * D_h**4) / (128 * mu_iso * Lk)
    
A_h = AssemblyHidraulico(conec, C_iso)
A_h_tilde = A_h.copy()
A_h_tilde[n_outlet, :] = 0.0
A_h_tilde[n_outlet, n_outlet] = 1.0

b_h = np.zeros(len(Xno))
b_h[n_inlet_0] = Q_in_0
b_h[n_inlet_175] = Q_in_175

pressures_iso = np.linalg.solve(A_h_tilde, b_h)
P_max_iso = np.max(pressures_iso)

print("=" * 80)
print("INICIALIZAÇÃO DO ACOMPLAMENTO HIDRÁULICO-TÉRMICO")
print("=" * 80)

print(f"Pressão Máxima Isotérmica Base calculada: {P_max_iso:.2f} Pa\n")

# -----------------------------------------------------------------------------
# EXERCÍCIO 1, PARTE 2
# -----------------------------------------------------------------------------

# Exercício teórico, deve ser colocado na apresentação.

# -----------------------------------------------------------------------------
# EXERCÍCIO 2, PARTE 2
# -----------------------------------------------------------------------------

print("=" * 80)
print("ANÁLISE DE INTERPOLAÇÃO 2D")
print("=" * 80)

Xno_base, conec = GeraGrafo(levels=3)
Xno = Xno_base * 0.001 
Xno[:, 1] = Xno[:, 1] + Ly/2  # Ajuste de offset de Y

# DEFINIÇÃO DA MALHA GROSSEIRA EXIGIDA (61 x 31)
Nx_coarse, Ny_coarse = 61, 31  
x_coarse = np.linspace(0, Lx, Nx_coarse)
y_coarse = np.linspace(0, Ly, Ny_coarse)
X_coarse, Y_coarse = np.meshgrid(x_coarse, y_coarse, indexing='ij')
pts_coarse = np.column_stack([X_coarse.ravel(), Y_coarse.ravel()])

métodos = ['linear', 'cubic', 'nearest']

# -------------------------------------------------------------------------
# CASO A: Base de Referência com Malha Exata (241 x 121)
# -------------------------------------------------------------------------
print(">> Processando Caso A: Solução Malha Exata (241 x 121)...")
Nx_ref, Ny_ref = 241, 121
d_max_ref = 0.0005

x_ref_coords = np.linspace(0, Lx, Nx_ref)
y_ref_coords = np.linspace(0, Ly, Ny_ref)
TB_ref = 10.0 + 20.0 * (x_ref_coords / Lx)
TT_ref = 10.0 + 20.0 * (x_ref_coords / Lx)

k_e_ref, k_n_ref = ObterCondutividadeFaces_ViaNos(Nx_ref, Ny_ref, Lx, Ly, Xno, conec, d_max_ref, k_0)
A_s_ref, b_s_ref = CriarSistemaSolidoCondutividadeVariavel(
    Nx_ref, Ny_ref, Lx, Ly, k_e_ref, k_n_ref, TL, TR, TB_ref, TT_ref, fonte_calor, R_incl, xincl, yincl, TC
)
T_solid_ref = spsolve(sparse.csr_matrix(A_s_ref), b_s_ref)

data_ref = T_solid_ref.reshape((Ny_ref, Nx_ref)).T
T_min_ref, T_max_ref = np.min(T_solid_ref), np.max(T_solid_ref)

figA_cont, axesA_cont = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
figA_cont.suptitle(f"Caso A (Malha Exata 241x121) - Contornos projetados na Grade {Nx_coarse}x{Ny_coarse}", fontweight='bold', fontsize=13)

figA_graf, axesA_graf = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
figA_graf.suptitle("Caso A (Malha Exata 241x121) - Temperaturas nos Nós Hidráulicos", fontweight='bold', fontsize=13)

for idx, método in enumerate(métodos):
    interp_ref = RegularGridInterpolator((x_ref_coords, y_ref_coords), data_ref, method=método, bounds_error=False, fill_value=None)
    
    # 1. Contorno
    temps_coarse = interp_ref(pts_coarse).reshape(Nx_coarse, Ny_coarse)
    temps_coarse = np.clip(temps_coarse, T_min_ref, T_max_ref)
    
    cp = axesA_cont[idx].contourf(x_coarse, y_coarse, temps_coarse.T, levels=40, cmap='jet')
    axesA_cont[idx].set_title(f"Interpolação: {método.capitalize()}", fontsize=11)
    axesA_cont[idx].set_xlabel("X (m)")
    axesA_cont[idx].set_ylabel("Y (m)")
    axesA_cont[idx].set_aspect('equal')
    
    # 2. Grafo Nodal
    T_nodes_hyd_ref = interp_ref(Xno)
    T_nodes_hyd_ref = np.clip(T_nodes_hyd_ref, T_min_ref, T_max_ref)
    
    sc = plot_nos_cromaticos_hidraulica(conec, Xno, T_nodes_hyd_ref, método, Nx_ref, Ny_ref, ax=axesA_graf[idx])
    axesA_graf[idx].set_title(f"Método: {método.capitalize()}", fontsize=11)
    axesA_graf[idx].set_aspect('equal')

cbarA_cont = figA_cont.colorbar(cp, ax=axesA_cont.ravel().tolist(), orientation='horizontal', shrink=0.4, aspect=40)
cbarA_cont.set_label('Temperatura (°C)', fontsize=11)

cbarA_graf = figA_graf.colorbar(sc, ax=axesA_graf.ravel().tolist(), orientation='horizontal', shrink=0.4, aspect=40)
cbarA_graf.set_label('Temperatura do Nó (°C)', fontsize=11)

plt.show()

# -------------------------------------------------------------------------
# CASO B: Solução Base com Malha Grosseira (61 x 31)
# -------------------------------------------------------------------------

Nx_61, Ny_61 = 61, 31
d_max_61 = 0.0005

x_61_coords = np.linspace(0, Lx, Nx_61)
y_61_coords = np.linspace(0, Ly, Ny_61)
TB_61 = 10.0 + 20.0 * (x_61_coords / Lx)
TT_61 = 10.0 + 20.0 * (x_61_coords / Lx)

k_e_61, k_n_61 = ObterCondutividadeFaces_ViaNos(Nx_61, Ny_61, Lx, Ly, Xno, conec, d_max_61, k_0)
A_s_61, b_s_61 = CriarSistemaSolidoCondutividadeVariavel(
    Nx_61, Ny_61, Lx, Ly, k_e_61, k_n_61, TL, TR, TB_61, TT_61, fonte_calor, R_incl, xincl, yincl, TC
)
T_solid_61 = spsolve(sparse.csr_matrix(A_s_61), b_s_61)

data_61 = T_solid_61.reshape((Ny_61, Nx_61)).T
T_min_61, T_max_61 = np.min(T_solid_61), np.max(T_solid_61)

figB_cont, axesB_cont = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
figB_cont.suptitle(f"Caso B (Malha Grosseira 61x31) - Contornos projetados na Grade {Nx_coarse}x{Ny_coarse}", fontweight='bold', fontsize=13)

figB_graf, axesB_graf = plt.subplots(1, 3, figsize=(15, 6), constrained_layout=True)
figB_graf.suptitle("Caso B (Malha Grosseira 61x31) - Temperaturas nos Nós Hidráulicos", fontweight='bold', fontsize=13)

for idx, método in enumerate(métodos):
    interp_61 = RegularGridInterpolator((x_61_coords, y_61_coords), data_61, method=método, bounds_error=False, fill_value=None)
    
    # 1. Contorno
    temps_coarse_61 = interp_61(pts_coarse).reshape(Nx_coarse, Ny_coarse)
    temps_coarse_61 = np.clip(temps_coarse_61, T_min_61, T_max_61)
    
    cp61 = axesB_cont[idx].contourf(x_coarse, y_coarse, temps_coarse_61.T, levels=40, cmap='jet')
    axesB_cont[idx].set_title(f"Interpolação: {método.capitalize()}", fontsize=11)
    axesB_cont[idx].set_xlabel("X (m)")
    axesB_cont[idx].set_ylabel("Y (m)")
    axesB_cont[idx].set_aspect('equal')
    
    # 2. Grafo Nodal
    T_nodes_hyd_61 = interp_61(Xno)
    T_nodes_hyd_61 = np.clip(T_nodes_hyd_61, T_min_61, T_max_61)
    
    sc61 = plot_nos_cromaticos_hidraulica(conec, Xno, T_nodes_hyd_61, método, Nx_61, Ny_61, ax=axesB_graf[idx])
    axesB_graf[idx].set_title(f"Método: {método.capitalize()}", fontsize=11)
    axesB_graf[idx].set_aspect('equal')

cbarB_cont = figB_cont.colorbar(cp61, ax=axesB_cont.ravel().tolist(), orientation='horizontal', shrink=0.4, aspect=40)
cbarB_cont.set_label('Temperatura (°C)', fontsize=11)

cbarB_graf = figB_graf.colorbar(sc61, ax=axesB_graf.ravel().tolist(), orientation='horizontal', shrink=0.4, aspect=40)
cbarB_graf.set_label('Temperatura do Nó (°C)', fontsize=11)

plt.show()

# Calcula a temperatura das arestas para o restante do seu código (Seção 3)
T_arestas_final, _ = calcular_temperatura_media_arestas(int(Nx_ref), int(Ny_ref), Lx, Ly, Xno, conec, T_solid_ref)
# -----------------------------------------------------------------------------
# EXERCÍCIO 3, PARTE 2
# -----------------------------------------------------------------------------

C_acoplado = atualiza_condutancias(conec, Xno, T_arestas_final, Area_canal)

A_h_mod = AssemblyHidraulico(conec, C_acoplado)
A_h_mod_tilde = A_h_mod.copy()
A_h_mod_tilde[n_outlet, :] = 0.0
A_h_mod_tilde[n_outlet, n_outlet] = 1.0

pressures_mod = np.linalg.solve(A_h_mod_tilde, b_h)
P_max_mod = np.max(pressures_mod)

print("-" * 80)
print(f"Pressão Máxima (Isotérmica):          {P_max_iso:.2f} Pa")
print(f"Pressão Máxima (Acoplada Corrigida):  {P_max_mod:.2f} Pa")
print(f"Variação da Perda de Carga Realizada: {((P_max_mod - P_max_iso)/P_max_iso)*100:.2f}%")
print("-" * 80)

# Renderização dos perfis de controle térmico do fluido
plot_arestas_cromaticas_hidraulics(conec, Xno, T_arestas_final, titulo="Distribuição Térmica Integrada nas Arestas da Rede (N=1000)")
print("\n")


# =========================================================================
# 3. EXERCÍCIOS DA PARTE 2 (SEÇÃO 4.3.3 - ANÁLISE PARAMÉTRICA DO SÓLIDO)
# =========================================================================
print("=" * 80)
print("3. EXERCÍCIOS DA PARTE 2 (SEÇÃO 4.3.3): ANÁLISE PARAMÉTRICA")
print("=" * 80)

# --- TAREFA 4.3.3.1: Estudo de Convergência de Malha e d_max ---
print("Iniciando varredura paramétrica de malhas e distâncias de corte:")
malhas = [(61, 31), (121, 61), (241, 121)]
d_max_valores = [0.00025, 0.0005, 0.0010]

print(f"{'Malha (Nx x Ny)':<18} | {'d_max (mm)':<12} | {'T_max (°C)':<12} | {'Tempo K (s)':<12} | {'Tempo Solv (s)':<12}")
print("-" * 80)

for (Nx, Ny) in malhas:
    x_coords = np.linspace(0, Lx, Nx)
    TB = 10.0 + 20.0 * (x_coords / Lx)
    TT = 10.0 + 20.0 * (x_coords / Lx)
    
    for d_max_val in d_max_valores:
        t_k_ini = time.perf_counter()
        
        # CORREÇÃO: Utilizar Nx, Ny e d_max_val dinâmicos do laço, e não os de referência (_ref)
        k_e, k_n = ObterCondutividadeFaces_ViaNos(Nx, Ny, Lx, Ly, Xno, conec, d_max_val, k_0)
        
        t_k_fim = time.perf_counter()
        
        t_s_ini = time.perf_counter()
        
        # CORREÇÃO: Passar k_e e k_n diretamente
        A_s, b_s = CriarSistemaSolidoCondutividadeVariavel(
            Nx, Ny, Lx, Ly, k_e, k_n, TL, TR, TB, TT, fonte_calor, R_incl, xincl, yincl, TC
        )
        T_solid_flat = spsolve(sparse.csr_matrix(A_s), b_s)
        t_s_fim = time.perf_counter()
        
        print(f"{f'{Nx} x {Ny}':<18} | {d_max_val*1000:<12.3f} | {np.max(T_solid_flat):<12.2f} | {t_k_fim - t_k_ini:<12.4f} | {t_s_fim - t_s_ini:<12.4f}")
        
print("-" * 80)

# --- ESPAÇO PARA A TAREFA 4.3.3.2: EFEITO DO TERMO FONTE PERTURBADO (S_p) ---
# TODO: Implementar a modificação do vetor de forças 'b_s' adicionando a contribuição
# das intensidades I_j das correntes associadas às ramificações (homogênea vs heterogênea).
print("\n[INFO] Espaço reservado para a implementação da Tarefa 4.3.3.2 (Termo Fonte Epigenético S_p).")
# =========================================================================
# IMPLEMENTAR AQUI O LAÇO PARA S_0 = [+/-1e5, +/-5e5, +/-1e6] E VARIAÇÃO DE I_j
# =========================================================================
import numpy as np
from scipy import sparse

# Importações dos seus arquivos base (CORRIGIDO: adicionado roda_exercicio_4)
from RedeHidraulica.functions import GeraGrafo, AssemblyVectorC, Assembly
from functionsM import BuildMatrizes_Eigen_Circular
from functionsHM import Monta_Matriz_Global_Acoplada, Resolve_Passo_Tempo, roda_exercicio_4, roda_exercicio_5

# ==============================================================================
# 1. PARÂMETROS FÍSICOS E ADIMENSIONAIS (Dados da Seção 5.2.5 do PDF)
# ==============================================================================
# Rede Hidráulica
mu = 5e-4              # Pa.s
complex_level = 3
mm_to_m = 0.001
H_k = 1000e-6          # Largura do canal (ex: 1000 um para o item 2)

# Membrana Elástica
R_fisico = 0.25e-2     # Raio físico: 0.25 cm em m
e_esp = 0.1e-3         # Espessura: 0.1 mm em m
sigma = 200.0          # Tensão: N/m
rho = 900.0            # Densidade: kg/m^3

# Parâmetros Adimensionais (Base)
Lx_ad, Ly_ad = 2.0, 2.0
R_ad = 1.0
sigma_ad = 1.0
rho_ad = 1.0
e_ad = 1.0
beta_ad = 0.1          # Fator de amortecimento

# Discretização e Simulação
Nx, Ny = 51, 51        # Malha base da membrana
h_ad = Lx_ad / (Nx - 1)
dt = 0.025             # Passo de tempo (adimensional)
t_max = 12.0           # Tempo final (adimensional)

p_inlet_dim = 10000.0  # Pressão de entrada (Pa) - Pode variar nos exercícios

# ==============================================================================
# 2. CONSTRUÇÃO DOS SISTEMAS ISOLADOS
# ==============================================================================
print("Construindo a Rede Hidráulica...")
Xno, conec = GeraGrafo(levels=complex_level)
Xno = Xno * mm_to_m

# CORREÇÃO AQUI: Passando H_k e mu para a função com assinatura atualizada
C = AssemblyVectorC(conec, Xno, H_k, mu) 

A_net = Assembly(conec, C) 
A_net_sparse = sparse.csr_matrix(A_net) # Converte para esparsa

# Definição dos nós de controle
n_inlet = 0
n_outlet = 5 # O PDF afirma que o Outlet principal da espinha é o nó 5
np_net = A_net_sparse.shape[0]

print("Construindo a Membrana Elástica...")
K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h_ad)
nm = Nx * Ny

# ==============================================================================
# 3. ACOPLAMENTO GLOBAL
# ==============================================================================
print("Montando a Matriz Global Acoplada...")
Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
    K_mem, M_mem, A_net_sparse, 
    Nx, Ny, Lx_ad, Ly_ad, R_ad, h_ad, 
    dt, beta_ad, sigma, rho, e_esp, R_fisico, 
    n_inlet, n_outlet
)

# Adimensionaliza a pressão de entrada
p_inlet_adim = p_inlet_dim / pref

# Arrays de estado iniciais (repouso)
w_n = np.zeros(nm)
v_n = np.zeros(nm)

tempos = np.arange(0, t_max + dt, dt)

# Listas preparadas para guardar o histórico (Úteis para os exercícios)
hist_p_out = []
hist_w_centro = []
hist_q_out = []
no_centro = int((Ny//2) * Nx + (Nx//2))

# ==============================================================================
# 4. LOOP TEMPORAL (TIME-STEPPING)
# ==============================================================================
print(f"Iniciando simulação transiente com {len(tempos)} passos...")

for t in tempos:
    
    # 4.1 Resolve o passo
    w_n, v_n, p_n = Resolve_Passo_Tempo(
        Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, np_net, n_inlet
    )
    
    # 4.2 Guarda os dados de interesse para os exercícios (re-dimensionalizando)
    p_out_fisico = p_n[n_outlet] * pref
    w_centro_fisico = w_n[no_centro] * (0.01 * R_fisico)
    
    # Vazão é q_outlet = h_ad^2 * \sum(v) -> Apenas para manter salvo!
    q_out_adim = h_ad**2 * np.sum(U_sparse[n_outlet, :].toarray().flatten() * v_n)
    
    hist_p_out.append(p_out_fisico)
    hist_w_centro.append(w_centro_fisico)
    hist_q_out.append(q_out_adim)

# ==============================================================================
# 5. DISPARO DO EXERCÍCIO 4
# ==============================================================================

#roda_exercicio_4(
#    Nx, Ny, h_ad, dt, mu, R_fisico, e_esp, sigma, rho, 
#    Lx_ad, Ly_ad, R_ad, sigma_ad, rho_ad, e_ad
#)

roda_exercicio_5(
    Nx, Ny, h_ad, dt, mu, R_fisico, e_esp, sigma, rho, 
    Lx_ad, Ly_ad, R_ad, sigma_ad, rho_ad, e_ad
)

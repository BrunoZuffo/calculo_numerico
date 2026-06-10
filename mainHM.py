import numpy as np
from scipy import sparse

# Importações dos seus arquivos base (Ajuste os nomes se necessário)
from functions import GeraGrafo, AssemblyVectorC, Assembly
from functionsM import BuildMatrizes_Eigen_Circular
from functionsHM import Monta_Matriz_Global_Acoplada, Resolve_Passo_Tempo

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
# Obs: O ideal aqui é ter certeza que AssemblyVectorC está usando o H_k se necessário,
# ou ajustá-la futuramente para aceitar H_k dinamicamente nos exercícios.
C = AssemblyVectorC(conec, Xno)
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
# 2. FUNÇÃO DO EXERCÍCIO 4 (Completa com análise de picos e plot)
# ==============================================================================
def roda_exercicio_4():
    from scipy.signal import find_peaks
    print("\n--- Iniciando Exercício 4 ---")
    
    # Parâmetros ESPECÍFICOS do Exercício 4
    H_k_ex4 = 2000e-6    # 2000 um
    beta_ad_ex4 = 0.0    # Sem amortecimento
    p_inlet_dim_ex4 = 0.0 # Pressão nula
    t_max_ex4 = 20.0
    
    # 1. Monta a Rede Hidráulica localmente
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001
    C = AssemblyVectorC(conec, Xno, H_k_ex4, mu) # Usa o H_k específico!
    A_net_sparse = sparse.csr_matrix(Assembly(conec, C))
    n_inlet, n_outlet = 0, 5 
    
    # 2. Monta a Membrana localmente
    K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h_ad)
    nm = Nx * Ny
    
    # 3. Acha o 3º Modo
    lambdas, Phi = eigsh(K_mem, k=5, M=M_mem, which='SM')
    idx = np.argsort(lambdas)
    lambdas, Phi = lambdas[idx], Phi[:, idx]
    modo_3 = Phi[:, 2]
    
    fator_Hz = (1.0 / (2.0 * np.pi * R_fisico)) * np.sqrt(sigma / (rho * e_esp))
    freq_analitica_Hz = np.sqrt(lambdas[2]) * fator_Hz
    print(f"-> Frequência Analítica do 3º Modo (Isolada): {freq_analitica_Hz:.2f} Hz")
    
    # 4. Acoplamento
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem, M_mem, A_net_sparse, Nx, Ny, Lx_ad, Ly_ad, R_ad, h_ad, 
        dt, beta_ad_ex4, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet
    )
    
    # 5. Loop no tempo
    p_inlet_adim = p_inlet_dim_ex4 / pref
    w_n = (modo_3 / np.max(np.abs(modo_3))) * 0.1 
    v_n = np.zeros(nm)
    
    tempos = np.arange(0, t_max_ex4 + dt, dt)
    hist_w_centro = []
    no_centro = int((Ny//2) * Nx + (Nx//2))
    
    for t in tempos:
        w_n, v_n, p_n = Resolve_Passo_Tempo(Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, A_net_sparse.shape[0], n_inlet)
        hist_w_centro.append(w_n[no_centro] * (0.01 * R_fisico))
        
    # 6. Processamento dos Picos e Cálculo da Frequência Simulada
    hist_w_array = np.array(hist_w_centro)
    picos, _ = find_peaks(hist_w_array)
    
    if len(picos) > 1:
        # Intervalo de tempo adimensional entre os picos encontrados
        tempos_picos = tempos[picos]
        periodos_ad = np.diff(tempos_picos)
        periodo_medio_ad = np.mean(periodos_ad)
        
        # Converte a frequência de oscilação para Hz
        freq_simulada_ad = 1.0 / periodo_medio_ad
        freq_simulada_Hz = freq_simulada_ad * fator_Hz
        
        # Erro percentual relativo
        erro_percentual = abs(freq_analitica_Hz - freq_simulada_Hz) / freq_analitica_Hz * 100
        
        print(f"-> Frequência Simulada (Acoplada): {freq_simulada_Hz:.2f} Hz")
        print(f"-> Diferença (Erro devido ao acoplamento): {erro_percentual:.2f}%")
    else:
        print("-> Erro: Não foram encontrados picos suficientes na oscilação.")
        freq_simulada_Hz = 0.0

    # 7. Plotagem dos Resultados do Exercício 4
    plt.figure(figsize=(10, 5))
    plt.plot(tempos, hist_w_centro, 'r-', linewidth=1.5, label='Deslocamento do Nó Central')
    
    if len(picos) > 1:
        plt.plot(tempos[picos], hist_w_array[picos], 'bo', markersize=6, label='Picos Identificados')
        
    plt.title(f"Exercício 4: Decaimento Livre Transiente (3º Modo)\n"
              f"Freq. Analítica: {freq_analitica_Hz:.2f} Hz | Freq. Simulada: {freq_simulada_Hz:.2f} Hz")
    plt.xlabel("Tempo Adimensional")
    plt.ylabel("Deslocamento Central (m)")
    plt.legend(loc="upper right")
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.show()

    

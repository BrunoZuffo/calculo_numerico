import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d

# --- Importações dos Módulos Anteriores ---
from RedeHidraulica.functions import GeraGrafo, Assembly
from AcoplamentoHidraulicoTermico.functionsHT import (
    ObterCondutividadeFaces_ViaNos,
    CriarSistemaSolidoCondutividadeVariavel,
    calcular_temperatura_media_arestas
)

# ==============================================================================
# LEIS CONSTITUTIVAS (SEÇÃO 6.2)
# ==============================================================================

# ---- Parametros fisicos  --

#Rede Hidráulica
H_k = 1000e-6          # Largura seção transversal quadrada
p_inlet = 5000.0      # Pressão no Nó de entrada
n_inlet = 0           # Nó de entrada
n_outlet = 5          # Nó de saída
complex_level = 3     # Nível de complexidade da rede

# Viscosidade dinâmica do fluido
def mu(T):
    return 0.001791 / (1 + 0.03368 * T + 0.000221 * (T**2))

#Placa termica
Lx = 0.03 
Ly = 0.015
k = 0.25               # Condutividade térmica
fonte = 5*pow(10,5)    # Fonte de calor
def h(Lx, Nx):
  return Lx / (Nx - 1)
TL= 10.0
TR = 30.0
def TB(x):
  return 10+20*(x/Lx)
def TT(x):
  return 10+20*(x/Lx)

#Temperatura circular:
xc = 0.0225
yc = 0.0075
R = 0.0025
TC = 35

#Membrana elástica:

R = 0.0025        # raio da membrana (0.25 cm)
e = 1.0e-4        # espessura (0.1 mm)
sigma = 200.0     # tensao membranal
rho = 900.0       # densidade
beta = 0.1        # Coeficiente de atrito viscoso (adimensional)


def Atualiza_Condutancias_GD(conec, Xno, L_k, T_arestas):
    """ Calcula a condutância dos canais com D_k e A_k rigorosos. """
    A_k = H_k * H_k  # Área da secção transversal: 500um x 500um
    D_k = np.sqrt(4 * A_k / np.pi)
    
    mu_k = mu(T_arestas)
    kappa_k = (np.pi * (D_k**4)) / (128 * mu_k)
    C_k = kappa_k / L_k
    return C_k


# ==============================================================================
# SEÇÃO 6.2 - BASE DO GEMEO DIGITAL E ESTADO NOMINAL
# ==============================================================================

def Setup_Base_GD():
    """ Configura a geometria estrita baseada na Secção 6.2. """
    Lx_placa, Ly_placa = 0.03, 0.015  
    Nx_term, Ny_term = 121, 61
    
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001
    Xno[:, 1] += 0.5 * Ly_placa 
    
    # Pré-cálculo dos comprimentos (L_k) das arestas para eficiência
    L_k = np.zeros(len(conec))
    for i, c in enumerate(conec):
        L_k[i] = np.linalg.norm(Xno[c[0]] - Xno[c[1]])
        
    return {
        'Xno': Xno, 'conec': conec, 'L_k': L_k,
        'Lx_placa': Lx_placa, 'Ly_placa': Ly_placa,
        'Nx_term': Nx_term, 'Ny_term': Ny_term
    }

def GD_Gera_Condutancias_Nominais(base):
    """ 
    Resolve a física Térmica apenas uma vez para gerar as condutâncias base.
    A temperatura não muda com as falhas hidráulicas.
    """
    Lx_p, Ly_p = base['Lx_placa'], base['Ly_placa']
    Nx_t, Ny_t = base['Nx_term'], base['Ny_term']
    Xno, conec = base['Xno'], base['conec']
    
    k0 = 0.25
    fonte = 5e5
    TC = 35.0
    
    d_max = 0.001
    k_e, k_n = ObterCondutividadeFaces_ViaNos(Nx_t, Ny_t, Lx_p, Ly_p, Xno, conec, d_max, k0)
    
    x_c = np.linspace(0, Lx_p, Nx_t)
    TB = 10.0 + 20.0 * (x_c / Lx_p)
    TT = 10.0 + 20.0 * (x_c / Lx_p)
    
    A_T, b_T = CriarSistemaSolidoCondutividadeVariavel(
        Nx_t, Ny_t, Lx_p, Ly_p, k_e, k_n, 
        10.0, 30.0, TB, TT, fonte, 0.0025, 0.0225, 0.0075, TC
    )
    T_flat = spsolve(sparse.csr_matrix(A_T), b_T)
    
    T_arestas, _ = calcular_temperatura_media_arestas(
        conec, Xno, T_flat, Nx_t, Ny_t, Lx_p, Ly_p, num_subintervalos=10
    )
    
    C_T_nom = Atualiza_Condutancias_GD(conec, Xno, base['L_k'], T_arestas)
    
    return C_T_nom


# ==============================================================================
# SEÇÃO 6.3.2 (EXERCÍCIO 1)
# ==============================================================================

def RandomFail(C_original, p_fail, f_obs):
    """ 
    Aplica obstruções estocásticas de forma vetorizada e rápida.
    """
    C_mod = C_original.copy()
    mask = np.random.rand(len(C_mod)) < p_fail
    C_mod[mask] /= f_obs
    return C_mod

def GD_Avalia_Vazao_Falhas(base, C_T_nom, p_inlet=5000.0, p_fail=0.05, f_obs=100.0):
    """ 
    Resolve o subsistema Hidráulico sob perturbações e devolve a vazão resultante.
    (Utilizado no Exercício 1 - Monte Carlo Estacionário).
    """
    # 1. Aplicação de falhas
    C_T = RandomFail(C_T_nom, p_fail, f_obs) if p_fail > 0.0 else C_T_nom.copy()
        
    # 2. Montagem do Sistema
    A_net = Assembly(base['conec'], C_T)
    A_net_tilde = A_net.copy()
    b_h = np.zeros(len(base['Xno']))
    
    # 3. Condições de Contorno (Pressão Fixa)
    # Nó 0: Pressão de entrada
    A_net_tilde[0, :] = 0.0
    A_net_tilde[0, 0] = 1.0
    b_h[0] = p_inlet
    
    # Nó 5: Pressão de saída (atmosférica)
    A_net_tilde[5, :] = 0.0
    A_net_tilde[5, 5] = 1.0
    b_h[5] = 0.0
    
    pressures = spsolve(sparse.csr_matrix(A_net_tilde), b_h)
    
    # 4. Cálculo da Vazão via Álgebra Linear
    q_inlet = A_net[0, :].dot(pressures)
    
    return q_inlet


# ==============================================================================
# SEÇÃO 6.3.2 (EXERCÍCIO 2)
# ==============================================================================




# ==============================================================================
# SEÇÃO 6.4.3 - APRENDIZADO DE MODELOS (INTERPOLAÇÃO E REGRESSÃO) (EXERCÍCIOS 1 E 2)
# ==============================================================================

def Interpola_Potencia_Linear(t_dados, P_dados):
    """ 
    Exercício 1: Interpolação Linear Local (Splines lineares).
    Retorna uma função avaliável que liga os dados por retas.
    """
    return interp1d(t_dados, P_dados, kind='linear')

def Interpola_Potencia_Cubica(t_dados, P_dados):
    """ 
    Exercício 2: Interpolação Cúbica Local (Splines cúbicos).
    Retorna uma função avaliável garantindo suavidade nas interfaces.
    """
    return interp1d(t_dados, P_dados, kind='cubic')
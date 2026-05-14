import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
from scipy.sparse.linalg import eigsh

# Importa a função construtora base
from functionsM import BuildMatrizes_Eigen_Circular, DecomposicaoModal

# INÍCIO ----------------------------------------------------------------------------------
# Parâmetros físicos gerais e fixos do problema (Seção 3.5 do PDF)

sigma = 200.0     # Tensão (N/m)
rho = 900.0       # Densidade (kg/m^3)
e_esp = 0.1e-3    # Espessura (m) -> 0.1 mm
R = 0.4e-2        # Raio (m) -> 0.4 cm

# Fator de conversão dimensional: f_k = (w_ad / 2pi) * (1/R) * sqrt(sigma / (rho * e_esp))
fator_Hz = (1.0 / (2.0 * np.pi * R)) * np.sqrt(sigma / (rho * e_esp))

# Parâmetros adimensionais para o motor algébrico
Lx_ad = 2.0
Ly_ad = 2.0
sigma_ad = 1.0
rho_ad = 1.0
e_ad = 1.0

# 3.5.1 Exercício 1 -------------------------------------------------------------------------------

N1_ex1 = 61
N2_ex1 = 61
h_ex1 = Lx_ad / (N1_ex1 - 1)

K_ex1, M_ex1 = BuildMatrizes_Eigen_Circular(N1_ex1, N2_ex1, sigma_ad, rho_ad, e_ad, h_ex1)

diag_K = K_ex1.diagonal()
mask_visual = diag_K.reshape((N2_ex1, N1_ex1))

plt.figure(figsize=(6,6))
plt.imshow(mask_visual > 1e5, extent=[0, Lx_ad, 0, Ly_ad], origin='lower', cmap='viridis')
plt.title(f"Validação do Exercício 1: Máscara da Membrana ({N1_ex1}x{N2_ex1})\n(Região Amarela = Pontos Fixos)")
plt.xlabel(r"$\hat{x}$")
plt.ylabel(r"$\hat{y}$")

plt.show() 

# 3.5.1 Exercício 2 -------------------------------------------------------------------------------

print("\n" + "="*115)
print(f"{'RESULTADOS DO EX2 (Convergência e Modos de Vibração)':^115}")
print("="*115)

casos_ex2 = [(21, 21), (41, 41), (61, 61), (81, 81), (101, 101)]
resultados_ex2 = []

Phi_final = None
omega_final = None
f_final = None

print(f"{'Malha':<10} | " + " | ".join([f"f{i+1:<7}" for i in range(10)]))
print("-" * 115)

for Nx, Ny in casos_ex2:
    h = Lx_ad / (Nx - 1)
    
    K, M = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h)
    
    # Extração de autovalores pelo método de Krylov
    lambdas, Phi = eigsh(K, k=10, M=M, which='SM')
    
    # Ordenação forçada (menor energia para maior)
    idx = np.argsort(lambdas)
    lambdas = lambdas[idx]
    Phi = Phi[:, idx]
    
    omega_ad = np.sqrt(np.abs(lambdas))
    f_Hz = omega_ad * fator_Hz
    
    resultados_ex2.append((Nx, f_Hz))
    
    linha = f"{Nx}x{Ny:<8} | " + " | ".join([f"{f:7.1f}" for f in f_Hz])
    print(linha)
    
    # Guarda as matrizes da malha mais refinada para a renderização 3D
    if Nx == 101:
        Phi_final = Phi
        omega_final = omega_ad
        f_final = f_Hz
        
print("=" * 115)

fig = plt.figure(figsize=(20, 8))
fig.canvas.manager.set_window_title('Exercício 2 - Modos de Vibração')

x = np.linspace(0, Lx_ad, 101)
y = np.linspace(0, Ly_ad, 101)
X, Y = np.meshgrid(x, y)

for i in range(10):
    ax = fig.add_subplot(2, 5, i+1, projection='3d')
    Z = Phi_final[:, i].reshape((101, 101))
    
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none', alpha=0.9)
    ax.set_title(f'Modo {i+1}\nf = {f_final[i]:.1f} Hz')
    
    z_limite = np.max(np.abs(Z))
    ax.set_zlim(-z_limite, z_limite)
    ax.axis('off')

plt.tight_layout()
plt.show()

# 3.5.1 Exercício 3 -------------------------------------------------------------------------------
# Tarefa puramente analítica. Será documentada em LaTeX no relatório final.
# Lembrete (Regra Rigorosa de Notação): a solução da oscilação livre deve ser descrita
# obrigatoriamente utilizando a notação sin, resultando na forma \sin(\omega_k t + \phi_k).


# 3.5.1 Exercício 4 -------------------------------------------------------------------------------
# (Placeholder) - Decomposição do termo forçante na base ortogonal de autovetores generalizados.

Nx = 101
Ny = 101

h = Lx_ad / (Nx - 1)

K, M = BuildMatrizes_Eigen_Circular(
    Nx, Ny,
    sigma_ad,
    rho_ad,
    e_ad,
    h
)

# modos naturais
lambdas, Phi = eigsh(K, k=10, M=M, which='SM')

idx = np.argsort(lambdas)

lambdas = lambdas[idx]
Phi = Phi[:, idx]

# -----------------------------------------------------------------
# construção do termo forçante espacial
# -----------------------------------------------------------------

x = np.linspace(0, Lx_ad, Nx)
y = np.linspace(0, Ly_ad, Ny)

X, Y = np.meshgrid(x, y)

Z_grid = (X - 0.5)**2 + (Y - 0.5)**2

Z = Z_grid.flatten()

# -----------------------------------------------------------------
# decomposição modal
# -----------------------------------------------------------------

alpha = DecomposicaoModal(Phi, M, Z)

print("\n" + "="*70)
print(f"{'DECOMPOSIÇÃO MODAL DO TERMO FORÇANTE':^70}")
print("="*70)

for i, a in enumerate(alpha):
    print(f"Modo {i+1:2d}: alpha = {a:.6e}")

print("="*70)


# 3.5.1 Exercício 5 -------------------------------------------------------------------------------
# (Placeholder) - Cálculo da energia elástica média e plotagem log-log da curva de ressonância.
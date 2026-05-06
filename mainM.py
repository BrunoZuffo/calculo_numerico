import numpy as np
from scipy.sparse.linalg import eigsh
import matplotlib
matplotlib.use("TkAgg")

from functionsM import BuildMatrizes_Eigen, PlotaModo

# ==============================================================================
# TESTE DA BASE DO CAPÍTULO 3 (SEÇÃO 3.4)
# Resolvendo os modos fundamentais de uma membrana quadrada de teste
# ==============================================================================

# Parâmetros Arbitrários (Malha e Geometria)
N1 = 41
N2 = 41
Lx = 1.0
Ly = 1.0
h = Lx / (N1 - 1)

# Propriedades Físicas (Arbitrárias para testar a consistência matricial)
sigma = 1.0  # Tensão
rho = 1.0    # Densidade
e = 1.0      # Espessura

# 1. Construção das Matrizes
print(f"Construindo matrizes K e M para malha de {N1}x{N2} nós...")
K, M = BuildMatrizes_Eigen(N1, N2, sigma, rho, e, h)

# 2. Resolução do Problema Generalizado de Autovalores
# Estamos interessados apenas nos modos de menor frequência (frequências fundamentais)
num_modos = 10
print(f"Extraindo os {num_modos} primeiros modos de vibração usando 'eigsh'...")

# eigsh resolve: K * Phi = lambda * M * Phi
# which='SM' -> Smallest in Magnitude (ignora os autovalores de penalidade "big_number")
lambdas, Phi = eigsh(K, k=num_modos, M=M, which='SM')

# As frequências naturais (rad/s) são a raiz quadrada dos autovalores
omegas = np.sqrt(np.abs(lambdas))

# 3. Impressão dos Resultados
print("\n" + "="*45)
print(f"{'FREQUÊNCIAS NATURAIS DE OSCILAÇÃO':^45}")
print("="*45)
print(f"{'Modo':>10} | {'Autovalor (lambda)':>20} | {'Frequência w (rad/s)':>20}")
print("-" * 55)

for i in range(num_modos):
    print(f"{i+1:10d} | {lambdas[i]:20.4f} | {omegas[i]:20.4f}")
print("=" * 55)

# 4. Plotagem Gráfica
# Visualizando a topologia dos 3 primeiros modos
for i in range(3):
    PlotaModo(N1, N2, Lx, Ly, Phi[:, i], i+1, omegas[i])
# --- SCRIPT DE TESTE PARA A SUA PARTE ---
# Importe as funções necessárias dos arquivos que você já tem
from functions import GeraGrafo
from functionsT import SolveSystemSparse
import numpy as np
from Testefuncao import TemperaturaMicrocanais

# 1. Gera a rede (mesmo código da hidráulica)
Xno, conec = GeraGrafo(levels=2) # Usa um nível menor para testar mais rápido
mm_to_m = 0.001
Xno = Xno * mm_to_m

# 2. Gera a temperatura (mesmo código da térmica)
Nx, Ny = 51, 26
Lx, Ly = 0.02, 0.01
h = Lx / (Nx - 1)
K = 0.25
TL, TR = 10, 30
fonte = 5.0e5
x_coords = np.linspace(0, Lx, Nx)
TB = 10 + 20 * (x_coords / Lx)
TT = 10 + 20 * (x_coords / Lx)

T_grid, _, _, _ = SolveSystemSparse(Nx, Ny, h, K, TL, TR, TB, TT, fonte)

# 3. CHAMA A SUA FUNÇÃO!
T_nos_canais = TemperaturaMicrocanais(Xno, conec, T_grid, Lx, Ly)

print("Temperaturas calculadas para os 5 primeiros canais:")
for i in range(5):
    print(f"Canal {i}: {T_nos_canais[i]:.2f} °C")
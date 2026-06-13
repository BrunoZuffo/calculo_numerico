import shutil
from pathlib import Path
from time import perf_counter

import numpy as np
import pandas as pd
from scipy import sparse

from functions import GeraGrafo, Assembly
from functionsM import BuildMatrizes_Eigen_Circular
from functionsHM import Monta_Matriz_Global_Acoplada, Resolve_Passo_Tempo


# ==============================================================================
# QUESTÃO 2 - CAPÍTULO 5
# SCRIPT 1: RODA TODAS AS SIMULAÇÕES E SALVA DADOS BRUTOS
# ==============================================================================

RAW_DIR = Path("resultados_ex2_cap5_raw")
LIMPAR_SAIDA_ANTIGA = True

if LIMPAR_SAIDA_ANTIGA and RAW_DIR.exists():
    shutil.rmtree(RAW_DIR)

RAW_DIR.mkdir(exist_ok=True)

SNAP_DIR = RAW_DIR / "snapshots_nominal"
SNAP_DIR.mkdir(exist_ok=True)


# ==============================================================================
# 1. PARÂMETROS DO ENUNCIADO
# ==============================================================================

# Rede hidráulica
mu = 5e-4
complex_level = 3
mm_to_m = 0.001

n_inlet = 0
n_outlet = 5

larguras_canal = [1000e-6, 1250e-6, 1500e-6, 1750e-6]
pressoes_inlet = [5e3, 1e4, 2e4]

# Membrana física
R_fisico = 0.25e-2      # 0.25 cm em m
e_esp = 0.1e-3          # 0.1 mm em m
sigma = 200.0           # N/m
rho = 900.0             # kg/m^3

# Membrana adimensional
Lx_ad = 2.0
Ly_ad = 2.0
R_ad = 1.0
sigma_ad = 1.0
rho_ad = 1.0
e_ad = 1.0
beta_ad = 0.1

# Malhas e passos pedidos
malhas = [(51, 51), (101, 101)]
dts = [0.00625, 0.0125, 0.025, 0.05]
t_final = 12.0

# Caso nominal usado para salvar perfis espaciais da membrana
H_nominal_um = 1000
p_nominal = 10000
Nx_nominal = 51
Ny_nominal = 51
dt_nominal = 0.025
tempos_snapshot = [0.0, 3.0, 6.0, 12.0]


# ==============================================================================
# 2. FUNÇÕES AUXILIARES
# ==============================================================================

def nome_caso(H_k, p_inlet, Nx, Ny, dt):
    H_um = int(round(H_k * 1e6))
    p_int = int(round(p_inlet))
    return f"H{H_um}um_p{p_int}Pa_N{Nx}x{Ny}_dt{dt}"


def calcula_condutancias_quadradas(conec, Xno, H_k, mu):
    """
    Calcula C_k = kappa_k / L_k para canais quadrados.

    A = H_k^2
    D_h = sqrt(4A/pi)
    kappa = pi D_h^4 / (128 mu)
    """

    conec = np.asarray(conec, dtype=int)
    nc = conec.shape[0]
    C = np.zeros(nc)

    area = H_k ** 2
    D_h = np.sqrt(4.0 * area / np.pi)
    kappa = np.pi * D_h**4 / (128.0 * mu)

    for k in range(nc):
        n1 = conec[k, 0]
        n2 = conec[k, 1]

        x1, y1 = Xno[n1]
        x2, y2 = Xno[n2]

        L_k = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        if L_k <= 0:
            raise ValueError(f"Aresta {k} tem comprimento nulo.")

        C[k] = kappa / L_k

    return C


# ==============================================================================
# 3. GERAÇÃO DA REDE HIDRÁULICA
# ==============================================================================

print("\nGerando rede hidráulica...")

Xno, conec = GeraGrafo(levels=complex_level)
Xno = Xno * mm_to_m
conec = np.asarray(conec, dtype=int)

print(f"Nós da rede: {Xno.shape[0]}")
print(f"Arestas da rede: {conec.shape[0]}")


# ==============================================================================
# 4. FUNÇÃO QUE RODA UM CASO
# ==============================================================================

def roda_caso(Nx, Ny, dt, p_inlet_dim, H_k, salvar_snapshots=False):
    h_ad = Lx_ad / (Nx - 1)
    nm = Nx * Ny
    no_centro = (Ny // 2) * Nx + (Nx // 2)

    # --------------------------------------------------------------------------
    # Rede hidráulica
    # --------------------------------------------------------------------------
    C = calcula_condutancias_quadradas(conec, Xno, H_k, mu)
    A_net = Assembly(conec, C)
    A_net_sparse = sparse.csr_matrix(A_net)
    np_net = A_net_sparse.shape[0]

    # --------------------------------------------------------------------------
    # Membrana
    # --------------------------------------------------------------------------
    K_mem, M_mem = BuildMatrizes_Eigen_Circular(
        Nx, Ny,
        sigma_ad, rho_ad, e_ad,
        h_ad
    )

    # --------------------------------------------------------------------------
    # Matriz global acoplada
    # --------------------------------------------------------------------------
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem,
        M_mem,
        A_net_sparse,
        Nx,
        Ny,
        Lx_ad,
        Ly_ad,
        R_ad,
        h_ad,
        dt,
        beta_ad,
        sigma,
        rho,
        e_esp,
        R_fisico,
        n_inlet,
        n_outlet
    )

    p_inlet_adim = p_inlet_dim / pref

    # Escalas físicas
    wref = 0.01 * R_fisico
    qref = (R_fisico ** 2) * vref
    dt_fisico = dt * R_fisico * np.sqrt((rho * e_esp) / sigma)

    # Estado inicial
    w_n = np.zeros(nm)
    v_n = np.zeros(nm)

    tempos = np.arange(0.0, t_final + 0.5 * dt, dt)

    U_out = U_sparse[n_outlet, :].toarray().ravel()

    hist_t = []
    hist_p_outlet = []
    hist_q_outlet = []
    hist_w_centro = []
    hist_w_max = []
    hist_volume = []
    hist_potencia = []

    volume_acumulado = 0.0
    snapshots_salvos = set()

    for t in tempos:
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob,
            M_mem,
            w_n,
            v_n,
            p_inlet_adim,
            dt,
            nm,
            np_net,
            n_inlet
        )

        p_outlet_fisico = p_n[n_outlet] * pref

        q_out_adim = h_ad**2 * np.sum(U_out * v_n)
        q_out_fisico = q_out_adim * qref

        w_centro_fisico = w_n[no_centro] * wref
        w_max_fisico = np.max(np.abs(w_n)) * wref

        volume_acumulado += q_out_fisico * dt_fisico
        potencia = p_inlet_dim * q_out_fisico

        hist_t.append(t)
        hist_p_outlet.append(p_outlet_fisico)
        hist_q_outlet.append(q_out_fisico)
        hist_w_centro.append(w_centro_fisico)
        hist_w_max.append(w_max_fisico)
        hist_volume.append(volume_acumulado)
        hist_potencia.append(potencia)

        # Salva perfil espacial da membrana só para o caso nominal.
        if salvar_snapshots:
            for ts in tempos_snapshot:
                if ts not in snapshots_salvos and abs(t - ts) <= 0.5 * dt:
                    wmap = w_n.reshape((Ny, Nx)) * wref
                    np.savetxt(
                        SNAP_DIR / f"w_t{ts:.2f}.csv",
                        wmap,
                        delimiter=","
                    )
                    snapshots_salvos.add(ts)

    historico = pd.DataFrame({
        "t_adimensional": hist_t,
        "p_outlet_Pa": hist_p_outlet,
        "q_outlet_m3s": hist_q_outlet,
        "w_centro_m": hist_w_centro,
        "w_max_m": hist_w_max,
        "volume_acumulado_m3": hist_volume,
        "potencia_W": hist_potencia,
    })

    resumo = {
        "H_um": int(round(H_k * 1e6)),
        "p_inlet_Pa": float(p_inlet_dim),
        "Nx": Nx,
        "Ny": Ny,
        "dt": dt,
        "n_passos": len(tempos),
        "nm_membrana": nm,
        "np_rede": np_net,
        "w_centro_max_m": float(np.max(np.abs(hist_w_centro))),
        "w_centro_final_m": float(hist_w_centro[-1]),
        "w_max_m": float(np.max(hist_w_max)),
        "p_outlet_max_Pa": float(np.max(hist_p_outlet)),
        "p_outlet_final_Pa": float(hist_p_outlet[-1]),
        "q_outlet_max_abs_m3s": float(np.max(np.abs(hist_q_outlet))),
        "q_outlet_final_m3s": float(hist_q_outlet[-1]),
        "volume_final_m3": float(hist_volume[-1]),
        "potencia_max_abs_W": float(np.max(np.abs(hist_potencia))),
        "potencia_final_W": float(hist_potencia[-1]),
    }

    return historico, resumo


# ==============================================================================
# 5. RODAR TODOS OS CASOS
# ==============================================================================

total_casos = len(larguras_canal) * len(pressoes_inlet) * len(malhas) * len(dts)
contador = 0
resumos = []

inicio_total = perf_counter()

print(f"\nTotal de casos: {total_casos}")

for H_k in larguras_canal:
    for p_inlet in pressoes_inlet:
        for Nx, Ny in malhas:
            for dt in dts:
                contador += 1

                caso = nome_caso(H_k, p_inlet, Nx, Ny, dt)

                print("\n" + "=" * 90)
                print(f"Rodando caso {contador}/{total_casos}: {caso}")
                print("=" * 90)

                salvar_snapshots = (
                    int(round(H_k * 1e6)) == H_nominal_um
                    and int(round(p_inlet)) == p_nominal
                    and Nx == Nx_nominal
                    and Ny == Ny_nominal
                    and abs(dt - dt_nominal) < 1e-12
                )

                inicio_caso = perf_counter()

                historico, resumo = roda_caso(
                    Nx=Nx,
                    Ny=Ny,
                    dt=dt,
                    p_inlet_dim=p_inlet,
                    H_k=H_k,
                    salvar_snapshots=salvar_snapshots
                )

                resumo["tempo_execucao_s"] = float(perf_counter() - inicio_caso)

                historico.to_csv(RAW_DIR / f"historico_{caso}.csv", index=False)
                resumos.append(resumo)

                print(f"Tempo do caso: {resumo['tempo_execucao_s']:.2f} s")

df_resumo = pd.DataFrame(resumos)
df_resumo = df_resumo.sort_values(["H_um", "p_inlet_Pa", "Nx", "dt"])
df_resumo.to_csv(RAW_DIR / "tabela_raw_ex2_cap5.csv", index=False)

tempo_total = perf_counter() - inicio_total

with open(RAW_DIR / "LEIA_ME_RAW.txt", "w", encoding="utf-8") as f:
    f.write(
        "Dados brutos da Questão 2 - Capítulo 5\n\n"
        f"Total de casos: {total_casos}\n"
        f"Tempo total de execução: {tempo_total:.2f} s\n\n"
        "Arquivos:\n"
        "- tabela_raw_ex2_cap5.csv: tabela completa com uma linha por caso.\n"
        "- historico_*.csv: evolução temporal de cada caso.\n"
        "- snapshots_nominal/w_t*.csv: perfis espaciais da membrana para o caso nominal.\n"
    )

print("\n" + "=" * 90)
print("SIMULAÇÃO PESADA FINALIZADA")
print("=" * 90)
print(f"Dados brutos salvos em: {RAW_DIR}")
print(f"Tabela principal: {RAW_DIR / 'tabela_raw_ex2_cap5.csv'}")
print(f"Tempo total: {tempo_total:.2f} s")
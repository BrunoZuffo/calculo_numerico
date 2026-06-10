import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy import sparse
from scipy.sparse.linalg import spsolve

from functions import (
    GeraGrafo,
    SolveNetwork,
    calc_vazao,
    calc_potencia
)

from functionsHT import (
    calcular_temperatura_media_arestas_campo,
    calcular_condutancias_por_temperatura,
    criar_interpolador_temperatura,
    CriarSistemaSolidoCondutividadeVariavel,
    plot_arestas_cromaticas_hidraulics
)


print("\n" + "="*70)
print("EXERCÍCIO 4 - CAPÍTULO 4: REDE HIDRÁULICA COM C(T)")
print("="*70)

# ---------------------------------------------------------------------
# Parâmetros da placa térmica do Capítulo 4
# ---------------------------------------------------------------------

Lx, Ly = 0.03, 0.015
k_0 = 0.25
fonte_calor = 5.0e5

TL, TR = 10.0, 30.0

R_incl = 0.0025
xincl, yincl = 0.02 + R_incl, 0.0075
TC = 35.0

# ---------------------------------------------------------------------
# Rede hidráulica do Capítulo 4
# ---------------------------------------------------------------------

Xno, conec = GeraGrafo(levels=3)
Xno = Xno * 0.001

# Centraliza a rede em y dentro da placa de altura Ly
Xno[:, 1] += 0.5 * Ly

# Pressão fixa no outlet.
# SolveNetwork usa índice começando em 1.
# Nó 5 em Python vira chave "6".
ps_hid = {
    "6": 0.0
}

# Vazões impostas:
# nó 0   -> chave "1"   -> 0.1 mL/s = 1e-7 m³/s
# nó 175 -> chave "176" -> 1.0 mL/s = 1e-6 m³/s
Qs_hid = {
    "1": 1.0e-7,
    "176": 1.0e-6
}

# ---------------------------------------------------------------------
# Casos: malha térmica grosseira/refinada e regras de quadratura
# ---------------------------------------------------------------------

lista_malhas_termicas = [
    (61, 31),
    (241, 121)
]

lista_metodos = [
    "ponto_medio",
    "trapezio"
]

lista_subintervalos = [
    1,
    10,
    100,
    1000
]

resultados_ex4 = []

ultima_T_arestas = None
ultimo_interpolador = None

for Nx, Ny in lista_malhas_termicas:
    print(f"\n--- Resolvendo campo térmico na malha {Nx} x {Ny} ---")

    x_coords = np.linspace(0.0, Lx, Nx)
    y_coords = np.linspace(0.0, Ly, Ny)

    TB = 10.0 + 20.0 * (x_coords / Lx)
    TT = 10.0 + 20.0 * (x_coords / Lx)

    # Condutividade térmica uniforme da placa
    k_e = np.ones((Nx - 1, Ny)) * k_0
    k_n = np.ones((Nx, Ny - 1)) * k_0

    # Montagem do sistema térmico
    t0 = time.perf_counter()

    A_T, b_T = CriarSistemaSolidoCondutividadeVariavel(
        Nx,
        Ny,
        Lx,
        Ly,
        k_e,
        k_n,
        TL,
        TR,
        TB,
        TT,
        fonte_calor,
        R_incl,
        xincl,
        yincl,
        TC
    )

    tempo_montagem_termica = time.perf_counter() - t0

    # Resolução térmica
    t0 = time.perf_counter()

    T_flat = spsolve(sparse.csr_matrix(A_T), b_T)

    tempo_solucao_termica = time.perf_counter() - t0

    # Interpolador T(x,y)
    interpolador_T, T_matriz = criar_interpolador_temperatura(
        T_flat,
        Nx,
        Ny,
        Lx,
        Ly,
        metodo="linear"
    )

    ultimo_interpolador = interpolador_T

    for metodo in lista_metodos:
        for nsub in lista_subintervalos:
            print(f"Malha {Nx}x{Ny} | Método: {metodo} | nsub = {nsub}")

            # 1. Temperatura média nas arestas
            T_arestas, tempo_quadratura = calcular_temperatura_media_arestas_campo(
                conec=conec,
                Xno=Xno,
                interpolador_T=interpolador_T,
                metodo=metodo,
                num_subintervalos=nsub
            )

            # 2. Condutâncias hidráulicas atualizadas por temperatura
            C_T = calcular_condutancias_por_temperatura(
                conec=conec,
                Xno=Xno,
                T_arestas=T_arestas,
                area_canal=2.5e-7
            )

            # 3. Resolver rede hidráulica com C(T)
            t0 = time.perf_counter()

            pressure_T = SolveNetwork(
                conec,
                C_T,
                ps=ps_hid,
                Qs=Qs_hid
            )

            tempo_hidraulico = time.perf_counter() - t0

            # 4. Vazões e potência
            Q_T = calc_vazao(conec, C_T, pressure_T)
            W_T = calc_potencia(conec, C_T, pressure_T)

            resultados_ex4.append({
                "malha_termica": f"{Nx}x{Ny}",
                "metodo_quadratura": metodo,
                "nsub": nsub,
                "T_aresta_min": np.min(T_arestas),
                "T_aresta_media": np.mean(T_arestas),
                "T_aresta_max": np.max(T_arestas),
                "pressao_max": np.max(pressure_T),
                "pressao_min": np.min(pressure_T),
                "potencia_total": W_T,
                "tempo_montagem_termica": tempo_montagem_termica,
                "tempo_solucao_termica": tempo_solucao_termica,
                "tempo_quadratura": tempo_quadratura,
                "tempo_hidraulico": tempo_hidraulico
            })

            ultima_T_arestas = T_arestas


# ---------------------------------------------------------------------
# Tabela final
# ---------------------------------------------------------------------

df_ex4 = pd.DataFrame(resultados_ex4)

print("\n" + "="*70)
print("RESULTADOS DO EXERCÍCIO 4")
print("="*70)
print(df_ex4)

df_ex4.to_csv("resultados_exercicio4_cap4.csv", index=False)
print("\nArquivo salvo: resultados_exercicio4_cap4.csv")


# ---------------------------------------------------------------------
# Gráfico 1: pressão máxima
# ---------------------------------------------------------------------

df_ex4["caso"] = (
    df_ex4["metodo_quadratura"]
    + " | N="
    + df_ex4["nsub"].astype(str)
)

plt.figure(figsize=(12, 5))

for malha in df_ex4["malha_termica"].unique():
    dados = df_ex4[df_ex4["malha_termica"] == malha]
    plt.plot(
        dados["caso"],
        dados["pressao_max"],
        marker="o",
        label=f"Malha {malha}"
    )

plt.title("Exercício 4 - Pressão máxima")
plt.xlabel("Regra de quadratura")
plt.ylabel("Pressão máxima")
plt.xticks(rotation=45, ha="right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("ex4_pressao_maxima.png", dpi=200)
plt.close()

print("Gráfico salvo: ex4_pressao_maxima.png")


# ---------------------------------------------------------------------
# Gráfico 2: potência total
# ---------------------------------------------------------------------

plt.figure(figsize=(12, 5))

for malha in df_ex4["malha_termica"].unique():
    dados = df_ex4[df_ex4["malha_termica"] == malha]
    plt.plot(
        dados["caso"],
        dados["potencia_total"],
        marker="o",
        label=f"Malha {malha}"
    )

plt.title("Exercício 4 - Potência total")
plt.xlabel("Regra de quadratura")
plt.ylabel("Potência total")
plt.xticks(rotation=45, ha="right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("ex4_potencia_total.png", dpi=200)
plt.close()

print("Gráfico salvo: ex4_potencia_total.png")


# ---------------------------------------------------------------------
# Visualização da rede colorida pela temperatura média das arestas
# ---------------------------------------------------------------------

T_nodes_plot = np.zeros(len(Xno))

for i in range(len(Xno)):
    T_nodes_plot[i] = ultimo_interpolador([Xno[i]])[0]

plot_arestas_cromaticas_hidraulics(
    conec,
    Xno,
    T_nodes=T_nodes_plot,
    T_arestas=ultima_T_arestas,
    titulo="Exercício 4 - Temperatura média nas arestas da rede hidráulica"
)

print("\nExercício 4 finalizado.")

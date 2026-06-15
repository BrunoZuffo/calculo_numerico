import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# ==============================================================================
# QUESTÃO 2 - CAPÍTULO 5
# SCRIPT 2: PÓS-PROCESSAMENTO
# Lê resultados_ex2_cap5_raw e gera arquivos finais em resultados_ex2_cap5_apresentacao
# ==============================================================================

RAW_DIR = Path("resultados_ex2_cap5_pre_processado")
OUT_DIR = Path("resultados_ex2_cap5_apresentacao")
LIMPAR_SAIDA_ANTIGA = True

if LIMPAR_SAIDA_ANTIGA and OUT_DIR.exists():
    shutil.rmtree(OUT_DIR)

OUT_DIR.mkdir(exist_ok=True)

# Aceita nomes diferentes para a tabela bruta, caso você tenha gerado com outro nome.
possiveis_tabelas = [
    RAW_DIR / "tabela_raw_ex2_cap5.csv",
    RAW_DIR / "tabela_completa_ex2_cap5.csv",
    RAW_DIR / "tabela_final_ex2_cap5.csv",
]

TABELA_RAW = None

for arq in possiveis_tabelas:
    if arq.exists():
        TABELA_RAW = arq
        break

if TABELA_RAW is None:
    raise FileNotFoundError(
        f"Não encontrei nenhuma tabela válida em {RAW_DIR}. "
        "Procure se existe algum arquivo .csv de tabela nessa pasta."
    )

print("Usando tabela bruta:")
print(TABELA_RAW)


# ==============================================================================
# 1. LER TABELA RAW E GERAR TABELA FINAL ÚNICA
# ==============================================================================

df = pd.read_csv(TABELA_RAW)
df = df.sort_values(["H_um", "p_inlet_Pa", "Nx", "dt"])

colunas_final = [
    "H_um",
    "p_inlet_Pa",
    "Nx",
    "Ny",
    "dt",
    "n_passos",
    "w_centro_max_m",
    "w_max_m",
    "p_outlet_max_Pa",
    "q_outlet_max_abs_m3s",
    "volume_final_m3",
    "potencia_max_abs_W",
    "tempo_execucao_s",
]

df_final = df[colunas_final].copy()
df_final.to_csv(OUT_DIR / "tabela_final_ex2_cap5.csv", index=False)

print("Tabela final gerada:")
print(OUT_DIR / "tabela_final_ex2_cap5.csv")


# ==============================================================================
# 2. CONFIGURAÇÕES DOS CASOS REPRESENTATIVOS
# ==============================================================================

H_NOM = 1000
P_NOM = 10000
DT_NOM = 0.025

DTS = sorted(df_final["dt"].unique())
DT_REF = min(DTS)

MALHAS = sorted(
    df_final[["Nx", "Ny"]].drop_duplicates().itertuples(index=False, name=None)
)


def nome_caso(H_um, p_inlet, Nx, Ny, dt):
    return f"H{int(H_um)}um_p{int(p_inlet)}Pa_N{int(Nx)}x{int(Ny)}_dt{dt}"


def ler_historico(H_um, p_inlet, Nx, Ny, dt):
    caso = nome_caso(H_um, p_inlet, Nx, Ny, dt)
    arquivo = RAW_DIR / f"historico_{caso}.csv"

    if not arquivo.exists():
        raise FileNotFoundError(f"Histórico não encontrado: {arquivo}")

    return pd.read_csv(arquivo)


def salvar_fig(nome):
    plt.tight_layout()
    plt.savefig(OUT_DIR / nome, dpi=220)
    plt.close()
    print("Figura gerada:", OUT_DIR / nome)


# ==============================================================================
# 3. GRÁFICO 01 - EVOLUÇÃO TEMPORAL COMPARANDO MALHAS
# ==============================================================================

variaveis = [
    ("w_max_m", "Deflexão máxima da membrana [m]"),
    ("p_outlet_Pa", "Pressão no outlet [Pa]"),
    ("q_outlet_m3s", "Vazão no outlet [m³/s]"),
    ("w_centro_m", "Deslocamento central [m]"),
    ("volume_acumulado_m3", "Volume acumulado [m³]"),
    ("potencia_W", "Potência hidráulica assinada [W]"),
]

fig, axs = plt.subplots(3, 2, figsize=(13, 12))
axs = axs.ravel()

for ax, (col, ylabel) in zip(axs, variaveis):
    for Nx, Ny in MALHAS:
        hist = ler_historico(H_NOM, P_NOM, Nx, Ny, DT_NOM)

        ax.plot(
            hist["t_adimensional"],
            hist[col],
            label=f"{Nx}x{Ny}"
        )

    ax.set_xlabel("Tempo adimensional")
    ax.set_ylabel(ylabel)
    ax.grid(True)

axs[0].legend(title="Malha")

fig.suptitle(
    f"Evolução temporal das grandezas pedidas | H={H_NOM} µm, "
    f"p_inlet={P_NOM} Pa, dt={DT_NOM}",
    fontsize=14
)

salvar_fig("01_evolucao_temporal_comparacao_malhas.png")


# ==============================================================================
# 4. GRÁFICO 02 - CONVERGÊNCIA TEMPORAL E MALHA
# ==============================================================================

df_nominal = df_final[
    (df_final["H_um"] == H_NOM) &
    (df_final["p_inlet_Pa"] == P_NOM)
].copy()

metricas_dt = [
    ("w_centro_max_m", "Deslocamento central máximo [m]"),
    ("p_outlet_max_Pa", "Pressão máxima no outlet [Pa]"),
    ("q_outlet_max_abs_m3s", "Vazão máxima no outlet [m³/s]"),
    ("volume_final_m3", "Volume final [m³]"),
]

fig, axs = plt.subplots(2, 2, figsize=(12, 9))
axs = axs.ravel()

for ax, (col, ylabel) in zip(axs, metricas_dt):
    for Nx, Ny in MALHAS:
        dados = df_nominal[
            (df_nominal["Nx"] == Nx) &
            (df_nominal["Ny"] == Ny)
        ].sort_values("dt")

        ax.plot(
            dados["dt"],
            dados[col],
            marker="o",
            label=f"{Nx}x{Ny}"
        )

    ax.set_xlabel("Passo de tempo adimensional dt")
    ax.set_ylabel(ylabel)
    ax.set_xticks(DTS)
    ax.set_xticklabels([str(dt) for dt in DTS])
    ax.grid(True)

axs[0].legend(title="Malha")
fig.suptitle("Convergência temporal e influência da malha", fontsize=14)

salvar_fig("02_convergencia_dt_malha.png")


# ==============================================================================
# 5. GRÁFICO 03 - COMPARAÇÃO ESPACIAL DIRETA ENTRE MALHAS
# ==============================================================================

df_malha = df_nominal[np.isclose(df_nominal["dt"], DT_REF)].copy()
df_malha["malha"] = df_malha["Nx"].astype(str) + "x" + df_malha["Ny"].astype(str)

metricas_malha = [
    ("w_centro_max_m", "Deslocamento central máximo [m]"),
    ("p_outlet_max_Pa", "Pressão máxima no outlet [Pa]"),
    ("q_outlet_max_abs_m3s", "Vazão máxima no outlet [m³/s]"),
    ("volume_final_m3", "Volume final [m³]"),
]

fig, axs = plt.subplots(2, 2, figsize=(12, 9))
axs = axs.ravel()

for ax, (col, ylabel) in zip(axs, metricas_malha):
    ax.bar(df_malha["malha"], df_malha[col])
    ax.set_xlabel("Malha")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y")

fig.suptitle(f"Comparação espacial direta entre malhas | dt={DT_REF}", fontsize=14)

salvar_fig("03_comparacao_direta_malhas.png")


# ==============================================================================
# 6. GRÁFICO 04 - SENSIBILIDADE À PRESSÃO DE ENTRADA COMPARANDO MALHAS
# ==============================================================================

metricas_pressao = [
    ("w_centro_max_m", "Deslocamento central máximo [m]"),
    ("p_outlet_max_Pa", "Pressão máxima no outlet [Pa]"),
    ("q_outlet_max_abs_m3s", "Vazão máxima no outlet [m³/s]"),
    ("potencia_max_abs_W", "Potência máxima absoluta [W]"),
]

fig, axs = plt.subplots(2, 2, figsize=(12, 9))
axs = axs.ravel()

for ax, (col, ylabel) in zip(axs, metricas_pressao):

    for Nx, Ny in MALHAS:
        df_pressao = df_final[
            (df_final["H_um"] == H_NOM) &
            (df_final["Nx"] == Nx) &
            (df_final["Ny"] == Ny) &
            (np.isclose(df_final["dt"], DT_NOM))
        ].copy()

        df_pressao = df_pressao.sort_values("p_inlet_Pa")

        ax.plot(
            df_pressao["p_inlet_Pa"],
            df_pressao[col],
            marker="o",
            label=f"{Nx}x{Ny}"
        )

    ax.set_xlabel("Pressão de entrada [Pa]")
    ax.set_ylabel(ylabel)
    ax.grid(True)

axs[0].legend(title="Malha")

fig.suptitle(
    f"Sensibilidade à pressão de entrada | H={H_NOM} µm, dt={DT_NOM}",
    fontsize=14
)

salvar_fig("04_sensibilidade_pressao.png")


# ==============================================================================
# 7. GRÁFICO 05 - SENSIBILIDADE À LARGURA DO CANAL COMPARANDO MALHAS
# ==============================================================================

metricas_largura = [
    ("w_centro_max_m", "Deslocamento central máximo [m]"),
    ("p_outlet_max_Pa", "Pressão máxima no outlet [Pa]"),
    ("q_outlet_max_abs_m3s", "Vazão máxima no outlet [m³/s]"),
    ("potencia_max_abs_W", "Potência máxima absoluta [W]"),
]

fig, axs = plt.subplots(2, 2, figsize=(12, 9))
axs = axs.ravel()

for ax, (col, ylabel) in zip(axs, metricas_largura):

    for Nx, Ny in MALHAS:
        df_largura = df_final[
            (df_final["p_inlet_Pa"] == P_NOM) &
            (df_final["Nx"] == Nx) &
            (df_final["Ny"] == Ny) &
            (np.isclose(df_final["dt"], DT_NOM))
        ].copy()

        df_largura = df_largura.sort_values("H_um")

        ax.plot(
            df_largura["H_um"],
            df_largura[col],
            marker="o",
            label=f"{Nx}x{Ny}"
        )

    ax.set_xlabel("Largura do canal [µm]")
    ax.set_ylabel(ylabel)
    ax.grid(True)

axs[0].legend(title="Malha")

fig.suptitle(
    f"Sensibilidade à largura do canal | p_inlet={P_NOM} Pa, dt={DT_NOM}",
    fontsize=14
)

salvar_fig("05_sensibilidade_largura.png")


# ==============================================================================
# 8. GRÁFICO 06 - PERFIL TRANSIENTE DA DEFLEXÃO DA MEMBRANA
# ==============================================================================

SNAP_DIR = RAW_DIR / "snapshots_nominal"


def extrai_tempo_snapshot(arq):
    # Exemplo: w_t3.00.csv -> 3.00
    return float(arq.stem.replace("w_t", ""))


snapshots = sorted(
    SNAP_DIR.glob("w_t*.csv"),
    key=extrai_tempo_snapshot
)

if snapshots:
    campos = []
    nomes = []

    for arq in snapshots:
        tempo = extrai_tempo_snapshot(arq)
        campos.append(np.loadtxt(arq, delimiter=","))
        nomes.append(f"t={tempo:.2f}")

    n = len(campos)

    fig, axs = plt.subplots(1, n, figsize=(4 * n, 4))

    if n == 1:
        axs = [axs]

    vmax = max(np.max(np.abs(campo)) for campo in campos)

    for ax, campo, nome in zip(axs, campos, nomes):
        im = ax.contourf(
            campo,
            levels=30,
            vmin=-vmax,
            vmax=vmax
        )
        ax.set_title(nome)
        ax.set_xlabel("i")
        ax.set_ylabel("j")

    fig.colorbar(im, ax=axs, shrink=0.85)
    fig.suptitle("Perfil transiente de deflexão da membrana | caso nominal", fontsize=14)

    salvar_fig("06_perfil_transiente_deflexao_membrana.png")

else:
    print("Aviso: não encontrei snapshots em", SNAP_DIR)
    print("O gráfico 06 não foi gerado.")


# ==============================================================================
# 9. ARQUIVO DE ORIENTAÇÃO PARA APRESENTAÇÃO
# ==============================================================================

texto = f"""
Questão 2 - Capítulo 5
Arquivos finais gerados em resultados_ex2_cap5_apresentacao

1) tabela_final_ex2_cap5.csv
Tabela única com todos os casos simulados.

2) 01_evolucao_temporal_comparacao_malhas.png
Mostra a evolução temporal das grandezas pedidas no enunciado:
deflexão máxima, pressão no outlet, vazão no outlet, deslocamento central,
volume acumulado e potência hidráulica assinada.

3) 02_convergencia_dt_malha.png
Mostra a influência do passo de tempo e da malha.

4) 03_comparacao_direta_malhas.png
Compara diretamente as malhas 51x51 e 101x101 usando o menor dt.

5) 04_sensibilidade_pressao.png
Mostra o efeito da pressão de entrada.
A potência mostrada nesse gráfico é a potência máxima absoluta.

6) 05_sensibilidade_largura.png
Mostra o efeito da largura dos canais.
A potência mostrada nesse gráfico é a potência máxima absoluta.

7) 06_perfil_transiente_deflexao_membrana.png
Mostra o perfil espacial da deflexão da membrana no caso nominal, em ordem temporal.

"""

with open(OUT_DIR / "LEIA_ME_APRESENTACAO.txt", "w", encoding="utf-8") as f:
    f.write(texto)


print("\n" + "=" * 90)
print("PÓS-PROCESSAMENTO FINALIZADO")
print("=" * 90)
print(f"Arquivos finais salvos em: {OUT_DIR}")
print("Use essa pasta na apresentação.")
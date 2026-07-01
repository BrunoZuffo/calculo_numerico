import os
import sys
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt

diretorio_atual = os.path.dirname(os.path.abspath(__file__))  
diretorio_pai = os.path.dirname(diretorio_atual)              

SUBPASTAS_PROJETO = [
    "RedeHidraulica",
    "PlacaTermica",
    "MembranaElastica",
    "AcoplamentoHidraulicoTermico",
    "AcoplamentoHidraulicoMecanico",
]

if diretorio_pai not in sys.path:
    sys.path.insert(0, diretorio_pai)

for _sub in SUBPASTAS_PROJETO:
    _caminho_sub = os.path.join(diretorio_pai, _sub)
    if os.path.isdir(_caminho_sub) and _caminho_sub not in sys.path:
        sys.path.insert(0, _caminho_sub)

from GemeoDigital.functionsGD import (
    MontaGemeoDigital, RodaSimulacaoBase, PlotaEstadoBaseGD,
    RodaExercicio1_FalhasHidraulicas, PlotaExercicio1_FalhasHidraulicas,
    RodaExercicio2_AnaliseDinamica, PlotaExercicio2_AnaliseDinamica,
    SimulaEnergiaDissipada, Interpola_Potencia_Linear, 
    Interpola_Potencia_Cubica
)

# ==============================================================================
# DEFINIÇÃO DOS PARÂMETROS INICIAIS DO SISTEMA
# ==============================================================================

params = {}

# ---- Conversoes e nivel de complexidade da rede ----
params['mm_to_m'] = 0.001
params['complex_level'] = 3          # Nivel de complexidade topologica: 3
params['n_inlet'] = 0                # No de entrada (Inlet 0)
params['n_outlet'] = 5               # No de saida (Outlet 5, espinha central)

# ---- Placa Termica ----
params['Lx'] = 0.03                          # Lx = 3 cm, em metros
params['Ly'] = 0.015                         # Ly = 1.5 cm, em metros
params['k0'] = 0.25                          # Condutividade da matriz solida: W/(K.m)
params['fonte_calor'] = 5.0e5                # Termo de fonte de calor: W/m^3
params['TL'] = 10.0                          # Temperatura no lado esquerdo: °C
params['TR'] = 30.0                          # Temperatura no lado direito: °C
params['R_incl'] = 0.25e-2                   # Raio da inclusao circular: 0.25 cm em m
params['xincl'] = 0.02 + params['R_incl']    # Centro em x: (2 + R) cm, em m
params['yincl'] = 0.75e-2                    # Centro em y: 0.75 cm, em m
params['TC'] = 35.0                          # Temperatura fixada na inclusao: °C

# Malha de discretizacao da placa termica e raio de influencia dos canais
# (compativel com as malhas usadas no acoplamento Hidraulico-Termico do Cap. 4)
params['Nx_t'] = 121
params['Ny_t'] = 61
params['d_max'] = 0.0005                     # Raio de proximidade (Equacao 4.3)

# ---- Rede Hidraulica ----
params['H_k'] = 1000e-6                      # Largura do canal (secao quadrada): 1000 um
params['p_inlet_dim'] = 5000.0               # Pressao prescrita no Inlet: 5000 Pa

# ---- Membrana Elastica ----
params['R_fisico'] = 0.25e-2                 # Raio fisico: 0.25 cm em m
params['e_esp'] = 0.1e-3                     # Espessura estrutural: 0.1 mm em m
params['sigma'] = 200.0                      # Tensao mecanica residual: N/m
params['rho'] = 900.0                        # Densidade de massa: kg/m^3
params['beta_ad'] = 0.1                      # Coeficiente de atrito viscoso (adimensional)

# Parametros adimensionais (base) e malha da membrana
params['Lx_ad'] = 2.0
params['Ly_ad'] = 2.0
params['R_ad'] = 1.0
params['sigma_ad'] = 1.0
params['rho_ad'] = 1.0
params['e_ad'] = 1.0
params['Nx_m'] = 51                          # Malha da membrana: 51 x 51 nos
params['Ny_m'] = 51

# Discretizacao temporal de verificacao
params['dt'] = 0.025
t_max_verificacao = 12.0

# ==============================================================================
# 2. MONTAGEM COMPLETA DO GEMEO DIGITAL
# ==============================================================================

print("=" * 80)
print("MONTAGEM DA BASE DO GEMEO DIGITAL (Capitulo 6)")
print("=" * 80)

print("Construindo Placa Termica, Rede Hidraulica e Membrana Elastica...")
GD = MontaGemeoDigital(params)


# =======================================================================================
# EXERCÍCIOS - CAPÍTULO 6
# =======================================================================================

# ==============================================================================
# SECAO 6.3.2 -  Investigando o comportamento do sistema via Monte Carlo
# ==============================================================================

# ========================================================================
# SECAO 6.3.2 - EXERCICIO 1
# Análise Estacionária de Falhas Hidráulicas
# ========================================================================

print("\n" + "=" * 80)
print("EXERCÍCIO 1 (Seção 6.3.2) - Análise Estacionária de Falhas Hidráulicas")
print("=" * 80)

# Limite crítico de vazão de entrada, conforme a "Pergunta de Projeto" do PDF
q_critico = 1.25e-5

N_convergencia = 8000
pO_convergencia = 0.6

# (c) Varredura no domínio de obstrução individual solicitado pelo PDF
pO_grid = np.arange(0.05, 0.65 + 1e-9, 0.05)
N_final = 3000   # tamanho amostral "convergido" usado em cada ponto da varredura

resultados_ex1 = RodaExercicio1_FalhasHidraulicas(
    GD,
    q_critico=q_critico,
    N_convergencia=N_convergencia,
    pO_convergencia=pO_convergencia,
    fObs_lista=(5, 10),
    pO_grid=pO_grid,
    N_final=N_final,
    p_outlet=0.0,
)

# ---- Salva os dois graficos do Exercicio 1 como imagem, na mesma pasta
#      deste script (diretorio_atual), alem de exibi-los na tela --------
caminho_fig_convergencia = os.path.join(diretorio_atual, "ex1_convergencia_monte_carlo.png")
caminho_fig_varredura = os.path.join(diretorio_atual, "ex1_prob_vs_pO.png")

PlotaExercicio1_FalhasHidraulicas(
    resultados_ex1,
    save_path_convergencia=caminho_fig_convergencia,
    save_path_varredura=caminho_fig_varredura,
)

print(f"\nFiguras do Exercicio 1 salvas em:")
print(f"  {caminho_fig_convergencia}")
print(f"  {caminho_fig_varredura}")

print(f"\nResumo da varredura (N = {resultados_ex1['N_final']} por ponto):")
for fObs in resultados_ex1['fObs_lista']:
    print(f"\n  f_obs = {fObs}")
    for pO, prob in zip(resultados_ex1['pO_grid'], resultados_ex1['prob_vs_pO'][fObs]):
        print(f"    p_O = {pO:.2f}  ->  Prob(q_inlet < q_crit) = {prob*100:.2f}%")

print("\n" + "=" * 80)
print("EXERCÍCIO 1 (Seção 6.3.2) CONCLUÍDO.")
print("=" * 80)

# ==============================================================================
# SEÇÃO 6.3.2 - EXERCÍCIO 2
# Análise Dinâmica do Gêmeo Digital Completo
# ==============================================================================

print("\n" + "=" * 80)
print("EXERCÍCIO 2 (Seção 6.3.2) - Análise Dinâmica do Gêmeo Digital")
print("=" * 80)

E_critico = 7.0
N_dinamica = 10  

resultados_ex2 = RodaExercicio2_AnaliseDinamica(
    GD,
    E_critico=E_critico,
    pO=pO_convergencia,          # mesmo cenário representativo do Exercício 1 (p_O = 0.6)
    fObs_lista=(5, 10),
    dt_lista=(0.05, 0.1),
    N=N_dinamica,
    t_max=4.0,
)

caminho_fig_ex2 = os.path.join(diretorio_atual, "ex2_analise_dinamica.png")
PlotaExercicio2_AnaliseDinamica(resultados_ex2, save_path=caminho_fig_ex2)

print(f"\nFigura do Exercício 2 salva em:\n  {caminho_fig_ex2}")

print(f"\nResumo do Exercício 2 (N = {resultados_ex2['N']} por combinação):")
for fObs in resultados_ex2['fObs_lista']:
    for dt in resultados_ex2['dt_lista']:
        prob = resultados_ex2['prob'][(fObs, dt)]
        print(f"  f_obs={fObs:<3d}  dt={dt:<5}  ->  "
              f"Prob(E < {E_critico}) = {prob * 100:.2f}%")

print("\n" + "=" * 80)
print("EXERCÍCIO 2 (Seção 6.3.2) CONCLUÍDO.")
print("=" * 80)


# ==============================================================================
# SECAO 6.4.3 Investigando o comportamento do sistema via aproximação de dados
# ==============================================================================

# ========================================================================
# SEÇÃO 6.4.3 - EXERCÍCIO 1
# Interpolação Linear Local
# SEÇÃO 6.4.3 - EXERCÍCIO 2
# Interpolação Cúbica Local
# ========================================================================

print("\n" + "-"*85)
print(" EXERCÍCIOS 1 E 2 (Seção 6.4.3): Interpolação linear e cúbca da potência")
print("-" * 85)

# 1. Definindo o vetor de tempo conforme exigido: t em [0, 4] com dt = 0.05
dt = 0.05
t_max_interp = 4.0

print("-> Simulando a trajetória nominal (sem falhas) do GD completo "
      f"para obter P(t) em [0, {t_max_interp}] com δt = {dt}...")
_E_nominal, P_dados, t_dados = SimulaEnergiaDissipada(
    GD, GD['C'], dt, t_max_interp, params['p_inlet_dim']
)
print(f"   Energia dissipada nominal: E = {_E_nominal:.4f}")

# 3. Gerando os modelos de interpolação
print("-> Construindo Splines Lineares (Ex 1)...")
interpolador_linear = Interpola_Potencia_Linear(t_dados, P_dados)

print("-> Construindo Splines Cúbicos (Ex 2)...")
interpolador_cubico = Interpola_Potencia_Cubica(t_dados, P_dados)

# 4. Avaliando os modelos em uma malha temporal bem fina para ver a "suavidade"
t_fino = np.linspace(t_dados[0], t_dados[-1], 500)
P_linear = interpolador_linear(t_fino)
P_cubico = interpolador_cubico(t_fino)


# ========================================================================
# 5. Renderização Gráfica Comparativa
# ========================================================================

# ------------------------------------------
# GRÁFICO 1: VISÃO GERAL (COMPLETO)
# ------------------------------------------
fig1, ax1 = plt.subplots(figsize=(10, 6))

# Pontos reais observados (Nós da malha)
ax1.plot(t_dados, P_dados, 'o', color='black', markersize=4, label='Dados Discretos Observados ($\\delta t = 0.05$)')

# Curvas de interpolação
ax1.plot(t_fino, P_linear, '--', color='dodgerblue', alpha=0.8, linewidth=1.5, label='Spline Linear (Ex 1)')
ax1.plot(t_fino, P_cubico, '-', color='crimson', alpha=0.8, linewidth=1.5, label='Spline Cúbico (Ex 2)')

# Formatação visual
ax1.set_title("Aproximação de Dados da Potência: Splines Lineares vs Cúbicos", fontweight='bold')
ax1.set_xlabel("Tempo adimensional ($t$)")
ax1.set_ylabel("Potência Instantânea $\\mathcal{P}(t)$")
ax1.legend(fontsize=11)
ax1.grid(True, linestyle='--', alpha=0.6)

fig1.tight_layout()
caminho_salvar_interp = os.path.join(diretorio_atual, "interpolacao_potencia_ex1_ex2.png")
fig1.savefig(caminho_salvar_interp, dpi=300)
print(f"✓ Gráfico de interpolação salvo em: {caminho_salvar_interp}")


# ------------------------------------------
# GRÁFICO 2: ZOOM NAS INTERFACES
# ------------------------------------------
fig2, ax2 = plt.subplots(figsize=(10, 6))

# Repetimos o plot na nova figura
ax2.plot(t_dados, P_dados, 'o', color='black', markersize=4, label='Dados Discretos Observados ($\\delta t = 0.05$)')
ax2.plot(t_fino, P_linear, '--', color='dodgerblue', alpha=0.8, linewidth=1.5, label='Spline Linear (Ex 1)')
ax2.plot(t_fino, P_cubico, '-', color='crimson', alpha=0.8, linewidth=1.5, label='Spline Cúbico (Ex 2)')

t_zoom_ini = t_dados[0] + 0.25 * (t_dados[-1] - t_dados[0])
t_zoom_fim = t_dados[0] + 0.50 * (t_dados[-1] - t_dados[0])
mascara_zoom = (t_fino >= t_zoom_ini) & (t_fino <= t_zoom_fim)
P_zoom = np.concatenate([P_cubico[mascara_zoom], P_linear[mascara_zoom]])
margem = 0.1 * (P_zoom.max() - P_zoom.min() + 1e-12)


ax2.set_xlim(1.4, 1.6)
ax2.set_ylim(0.0, 0.5)




# Formatação visual do zoom
ax2.set_title(f"Zoom ({1.4} a {1.6}s) - Analisando Suavidade nas Interfaces", fontweight='bold')
ax2.set_xlabel("Tempo adimensional ($t$)")
ax2.set_ylabel("Potência Instantânea $\\mathcal{P}(t)$")
ax2.legend(fontsize=11)
ax2.grid(True, linestyle='--', alpha=0.6)

fig2.tight_layout()
caminho_salvar_zoom = os.path.join(diretorio_atual, "interpolacao_potencia_zoom.png")
fig2.savefig(caminho_salvar_zoom, dpi=300)
print(f"✓ Gráfico com zoom salvo em: {caminho_salvar_zoom}")


plt.show()

print("\n" + "=" * 80)
print("EXERCÍCIOS 1 e 2 (Seção 6.4.3) CONCLUÍDO.")
print("=" * 80)

# ========================================================================
# # SEÇÃO 6.4.3 - EXERCÍCIO 3
# Regressão Polinomial Global
# ========================================================================

print("\n" + "-"*85)
print(" EXERCÍCIO 3 (6.4.3): Regressão Polinomial Global (grau m ∈ [3, 15])")
print("-" * 85)

graus_ex3 = list(range(3, 16))

print(f"-> Ajustando regressão polinomial global para graus m = {graus_ex3[0]} a {graus_ex3[-1]}...")
resultados_ex3 = Regressao_Polinomial_Global(t_dados, P_dados, graus=graus_ex3)

# ---- Tabela de erros L2 ----
print(f"\n{'Grau m':>8}  {'Erro L2 (Eq. 6.14)':>20}")
print("-" * 32)
for m in graus_ex3:
    print(f"{m:>8d}  {resultados_ex3['erros_L2'][m]:>20.6e}")

# ---- Figura 1: curvas de aproximação sobrepostas aos dados ----
t_fino_ex3 = np.linspace(t_dados[0], t_dados[-1], 600)

fig3, ax3 = plt.subplots(figsize=(12, 6))
ax3.plot(t_dados, P_dados, 'ko', markersize=3, label='Dados discretos ($\\delta t = 0.05$)', zorder=5)

cmap_ex3 = plt.cm.rainbow(np.linspace(0, 1, len(graus_ex3)))
for cor, m in zip(cmap_ex3, graus_ex3):
    coefs = resultados_ex3['coeficientes'][m]
    ax3.plot(t_fino_ex3, np.polyval(coefs, t_fino_ex3),
             color=cor, linewidth=1.2, alpha=0.85, label=f'm = {m}')

ax3.set_title("Regressão Polinomial Global — Graus m ∈ [3, 15]", fontweight='bold')
ax3.set_xlabel("Tempo adimensional ($t$)")
ax3.set_ylabel("Potência Instantânea $\\mathcal{P}(t)$")
ax3.legend(fontsize=8, ncol=3, loc='upper right')
ax3.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

caminho_fig_ex3 = os.path.join(diretorio_atual, "ex3_regressao_polinomial.png")
plt.savefig(caminho_fig_ex3, dpi=300)
print(f"\n✓ Figura das curvas de regressão salva em: {caminho_fig_ex3}")

# ---- Figura 2: erro L2 vs grau ----
fig3b, ax3b = plt.subplots(figsize=(8, 5))
erros_lista = [resultados_ex3['erros_L2'][m] for m in graus_ex3]
ax3b.semilogy(graus_ex3, erros_lista, 'o-', color='steelblue', linewidth=2, markersize=6)
ax3b.set_title("Erro L2 vs. Grau do Polinômio (Regressão Global)", fontweight='bold')
ax3b.set_xlabel("Grau $m$")
ax3b.set_ylabel("Erro $e$ (Eq. 6.14) — escala log")
ax3b.set_xticks(graus_ex3)
ax3b.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()

caminho_fig_ex3b = os.path.join(diretorio_atual, "ex3_erro_L2_vs_grau.png")
plt.savefig(caminho_fig_ex3b, dpi=300)
print(f"✓ Figura do erro L2 salva em: {caminho_fig_ex3b}")

plt.show()

# ========================================================================
# SEÇÃO 6.4.3 - EXERCÍCIO 4
# Análise de Sensibilidade a Ruídos Estocásticos
# ========================================================================

print("\n" + "-"*85)
print(" EXERCÍCIO 4 (6.4.3): Sensibilidade a Ruídos em p_inlet(t)")
print("-" * 85)

N_realizacoes_ex4 = 10   # número de trajetórias ruidosas
graus_ex4 = list(range(3, 16))

print(f"-> Simulando {N_realizacoes_ex4} trajetórias com p_inlet ruidosa "
      f"[U(-0.15, 0.15)] e ajustando regressão polinomial (m ∈ [{graus_ex4[0]}, {graus_ex4[-1]}])...")

resultados_ex4 = Analise_Ruido_Estocastico(
    GD,
    dt=dt,
    t_max=t_max_interp,
    p_inlet_dim=params['p_inlet_dim'],
    graus=graus_ex4,
    N_realizacoes=N_realizacoes_ex4,
    seed=42,
)

tempos_ex4 = resultados_ex4['tempos']

# ---- Tabela comparativa: erro médio com e sem ruído ----
print(f"\n{'Grau m':>8}  {'Erro L2 (sem ruído)':>22}  {'Erro L2 médio (com ruído)':>26}")
print("-" * 62)
for m in graus_ex4:
    e_limpo = resultados_ex3['erros_L2'][m]
    e_ruidoso = resultados_ex4['erros_L2_medio'][m]
    print(f"{m:>8d}  {e_limpo:>22.6e}  {e_ruidoso:>26.6e}")

# ---- Figura 1: trajetórias ruidosas + regressão de cada realização ----
t_fino_ex4 = np.linspace(tempos_ex4[0], tempos_ex4[-1], 600)

fig4, axes4 = plt.subplots(1, 2, figsize=(14, 5))

# Painel esquerdo — sinal ruidoso vs. dado limpo + média
ax_esq = axes4[0]
for k, P_r in enumerate(resultados_ex4['P_realizacoes']):
    ax_esq.plot(tempos_ex4, P_r, color='gray', linewidth=0.6, alpha=0.5,
                label='Realização ruidosa' if k == 0 else None)
ax_esq.plot(tempos_ex4, resultados_ex4['P_ruidoso_medio'],
            color='darkorange', linewidth=2.0, label='Média das realizações')
ax_esq.plot(t_dados, P_dados, 'k--', linewidth=1.5, label='Dados limpos (nominal)')
ax_esq.set_title("Sinal de Potência com Ruído Estocástico", fontweight='bold')
ax_esq.set_xlabel("Tempo adimensional ($t$)")
ax_esq.set_ylabel("Potência Instantânea $\\mathcal{P}(t)$")
ax_esq.legend(fontsize=9)
ax_esq.grid(True, linestyle='--', alpha=0.5)

# Painel direito — erro L2 médio (ruidoso) vs. grau, comparado ao limpo
ax_dir = axes4[1]
erros_limpos = [resultados_ex3['erros_L2'][m] for m in graus_ex4]
erros_ruidosos_med = [resultados_ex4['erros_L2_medio'][m] for m in graus_ex4]
ax_dir.semilogy(graus_ex4, erros_limpos, 's--', color='steelblue',
                linewidth=2, markersize=6, label='Sem ruído')
ax_dir.semilogy(graus_ex4, erros_ruidosos_med, 'o-', color='crimson',
                linewidth=2, markersize=6, label='Com ruído (média)')
ax_dir.set_title("Erro L2 vs. Grau: Sensibilidade ao Ruído", fontweight='bold')
ax_dir.set_xlabel("Grau $m$")
ax_dir.set_ylabel("Erro $e$ (Eq. 6.14) — escala log")
ax_dir.set_xticks(graus_ex4)
ax_dir.legend(fontsize=10)
ax_dir.grid(True, linestyle='--', alpha=0.5)

plt.tight_layout()
caminho_fig_ex4 = os.path.join(diretorio_atual, "ex4_sensibilidade_ruido.png")
plt.savefig(caminho_fig_ex4, dpi=300)
print(f"\n✓ Figura do Exercício 4 salva em: {caminho_fig_ex4}")
plt.show()

print("\n" + "=" * 80)
print("EXERCÍCIO 4 (Seção 6.4.3) CONCLUÍDO.")
print("=" * 80)


# ==============================================================================
# SEÇÃO 6.4.5 Investigando o comportamento do sistema via diferenciação numérica
# ==============================================================================

# ========================================================================
# SEÇÃO 6.4.5 - EXERCÍCIO 1
# Análise Numérica de Sensibilidade
# ========================================================================
print("\n" + "-"*85)
print(" EXERCÍCIO 1 (Seção 6.4.5): Análise Numérica de Sensibilidade")
print("-" * 85)

resultados_ex645 = RodaExercicio645_Sensibilidade(
    params,
    t_max=4.0,
    dt=0.05
)

PlotaExercicio645_Sensibilidade(
    resultados_ex645,
    save_dir=diretorio_atual
)

print("\nResumo numérico - Sensibilidade em relação a TC:")
for i, TC in enumerate(resultados_ex645["TC"]["valores"]):
    print(
        f"  TC={TC:8.2f} °C | "
        f"E={resultados_ex645['TC']['E'][i]:.6e} | "
        f"qin(tf)={resultados_ex645['TC']['qin_tf'][i]:.6e} | "
        f"V(tf)={resultados_ex645['TC']['V_tf'][i]:.6e}"
    )

print("\nResumo numérico - Sensibilidade em relação a H:")
for i, H in enumerate(resultados_ex645["H"]["valores"]):
    print(
        f"  H={H:8.2f} µm | "
        f"E={resultados_ex645['H']['E'][i]:.6e} | "
        f"dE/dH centered={resultados_ex645['H']['dE_centered'][i]:.6e} | "
        f"dE/dH direto={resultados_ex645['H']['dE_direct'][i]:.6e}"
    )

print("\n" + "=" * 80)
print("EXERCÍCIO 1 (Seção 6.4.5) CONCLUÍDO.")
print("=" * 80)

# ========================================================================
# SEÇÃO 6.4.5 - EXERCÍCIO 2
# Localização de Raízes via Newton-Raphson
# ========================================================================

print("\n" + "-"*85)
print(" EXERCÍCIO 2 (Seção 6.4.5): Localização de Raízes via Newton-Raphson")
print("-" * 85)

# parte do antero vai aqui

print("\n" + "=" * 80)
print("EXERCÍCIO 2 (Seção 6.4.5) CONCLUÍDO.")
print("=" * 80)

# =======================================================================================
# CONCLUSÃO
# =======================================================================================

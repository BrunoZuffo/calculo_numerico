"""
mainGD.py
----------------------------------------------------------------------
CAPITULO 6 - O GEMEO DIGITAL (GD)

Este script monta a BASE COMPLETA do Gemeo Digital, integrando os tres
subsistemas fisicos desenvolvidos ao longo da disciplina:

    1) Placa Termica      (Capitulo 2 / Capitulo 4 - parte 2)
    2) Rede Hidraulica     (Capitulo 1 / Capitulo 4 - parte 1)
    3) Membrana Elastica   (Capitulo 3)

acoplados monoliticamente conforme a abordagem do Capitulo 5, com as
especificacoes nominais definidas na Secao 6.2 do PDF da disciplina.

IMPORTANTE: este script resolve apenas a MONTAGEM e VERIFICACAO da
base do sistema (Secoes 6.1 e 6.2). Os exercicios de investigacao das
Secoes 6.3 (Predicoes sob Incerteza - Amostragem Monte Carlo) e 6.4
(Aprendizado de Modelos - Interpolacao, Minimos Quadrados e Zeros de
Funcao) devem ser implementados em sequencia, reaproveitando o
dicionario `GD` e a funcao `RodaSimulacaoBase` montados aqui. Pontos de
entrada sugeridos estao marcados ao final do arquivo.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("TkAgg")

# Garantindo que TODAS as pastas-irmas (RedeHidraulica, PlacaTermica,
# MembranaElastica, AcoplamentoHidraulicoTermico, AcoplamentoHidraulicoMecanico)
# fiquem diretamente visiveis no sys.path. Isso e necessario porque
# functionsHM.py (e possivelmente outros modulos) fazem imports
# "achatados" internamente (ex.: "from functions import ...",
# "from functionsM import ..."), que so resolvem se a pasta de cada
# modulo estiver diretamente no path - nao basta adicionar apenas a
# pasta-mae calculo_numerico.
diretorio_atual = os.path.dirname(os.path.abspath(__file__))   # .../GemeoDigital
diretorio_pai = os.path.dirname(diretorio_atual)                # .../calculo_numerico

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
)

# ==============================================================================
# 1. ESPECIFICACOES NOMINAIS DOS SUBSISTEMAS (Secao 6.2 do PDF)
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

#   print(f"\nRede hidraulica gerada com {GD['Xno'].shape[0]} nos e {GD['conec'].shape[0]} canos.")
#   print(f"Temperatura minima/maxima na placa: {GD['T_solid_flat'].min():.2f} °C / {GD['T_solid_flat'].max():.2f} °C")
#   print(f"Viscosidade minima/maxima nos canais: {GD['mu_arestas'].min():.4e} / {GD['mu_arestas'].max():.4e} Pa.s")
#   print(f"Condutancia minima/maxima dos canais: {GD['C'].min():.4e} / {GD['C'].max():.4e}")
#   print(f"Pressao de referencia (pref): {GD['pref']:.4f} Pa")
#   print(f"Velocidade de referencia (vref): {GD['vref']:.4e}")
#   
#   # ==============================================================================
#   # 3. SIMULACAO TRANSIENTE DE VERIFICACAO
#   # ==============================================================================
#   print("\nExecutando simulacao transiente de verificacao da base do GD...")
#   historico = RodaSimulacaoBase(GD, t_max=t_max_verificacao, p_inlet_dim=params['p_inlet_dim'])
#   
#   print(f"Pressao final no Outlet: {historico['p_out'][-1]:.4f} Pa")
#   print(f"Vazao adimensional final no Outlet: {historico['q_out'][-1]:.4e}")
#   print(f"Deslocamento final no centro da membrana: {historico['w_centro'][-1]*1e6:.4f} µm")
#   
#   # ==============================================================================
#   # 4. VISUALIZACAO DE VERIFICACAO
#   # ==============================================================================
#   PlotaEstadoBaseGD(GD, historico)
#   
#   print("\n" + "=" * 80)
#   print("BASE DO GEMEO DIGITAL MONTADA E VERIFICADA COM SUCESSO.")
#   print("Pronta para os exercicios das Secoes 6.3 (Monte Carlo) e")
#   print("6.4 (Interpolacao / Regressao / Zeros de Funcao).")
#   print("=" * 80)

# ==============================================================================
# 5. EXERCÍCIO 1 (Seção 6.3.2) - ANÁLISE ESTACIONÁRIA DE FALHAS HIDRÁULICAS
# ==============================================================================
#
# Subsistema ISOLADO (apenas Rede Hidráulica + Placa Térmica, em regime
# permanente, sem a membrana). Pressão de descarga nula no Outlet
# (p_outlet = 0). Pergunta de projeto: qual a probabilidade de que,
# após o desgaste operacional (obstrução estocástica dos microcanais),
# a vazão total de entrada q_inlet (nó 0) fique abaixo do limite
# crítico 1.25e-5 ?
# ==============================================================================
print("\n" + "=" * 80)
print("EXERCÍCIO 1 (Seção 6.3.2) - Análise Estacionária de Falhas Hidráulicas")
print("=" * 80)

# Limite crítico de vazão de entrada, conforme a "Pergunta de Projeto" do PDF
q_critico = 1.25e-5

# (b) Estudo de convergência: progressão de Prob(N) para um cenário
#     representativo de obstrução (p_O = 0.6, mesmo valor do exemplo
#     ilustrativo do PDF), contrastando f_obs = 5 e f_obs = 10.
#
#     N_convergencia = 8000 (em vez de 5000): com margem_alvo = 0.01
#     (1 ponto percentual, padrao da funcao DeterminaTamanhoAmostral-
#     Estabilizacao), o N necessario para phat ~ 0.835 (f_obs=5) fica
#     em torno de 5300 realizacoes - ou seja, com N_convergencia=5000
#     o criterio NUNCA era de fato atingido dentro da amostra, e o
#     N_min reportado ficava artificialmente preso no limite superior
#     (N_min = N_convergencia), nao representando uma estabilizacao
#     real. Com N_convergencia = 8000 ha margem suficiente para os
#     dois fatores de severidade (f_obs=5 e f_obs=10) estabilizarem
#     dentro do intervalo simulado, na ordem de grandeza de "algumas
#     milhares de realizacoes" mencionada no PDF.
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
# PRÓXIMOS PASSOS: DEMAIS EXERCÍCIOS DE INVESTIGAÇÃO
#
#   Seção 6.3.2, Item 2 - Análise Dinâmica do Gêmeo Digital Completo
#       (acoplamento forte multifísico, malha 51x51, indicador
#       energético E = integral de P(t) dt, malha temporal dt=0.05/0.1)
#
#   Seção 6.4.3 / 6.4.5 - Investigando o comportamento via aproximação
#       de dados e diferenciação numérica
#       - Reaproveitar GD['Aglob'], RodaSimulacaoBase(...) variando H_k
#         para gerar a família de curvas p(t; H) e ajustar/derivar A(H)
# ==============================================================================
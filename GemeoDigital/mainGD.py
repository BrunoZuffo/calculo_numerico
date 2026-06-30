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
N_dinamica = 2000  

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

# 5. Renderização Gráfica Comparativa
fig, ax = plt.subplots(figsize=(10, 6))

# Pontos reais observados (Nós da malha)
ax.plot(t_dados, P_dados, 'o', color='black', markersize=4, label='Dados Discretos Observados ($\\delta t = 0.05$)')

# Curvas de interpolação
ax.plot(t_fino, P_linear, '--', color='dodgerblue', alpha=0.8, linewidth=1.5, label='Spline Linear (Ex 1)')
ax.plot(t_fino, P_cubico, '-', color='crimson', alpha=0.8, linewidth=1.5, label='Spline Cúbico (Ex 2)')

# Formatação visual
ax.set_title("Aproximação de Dados da Potência: Splines Lineares vs Cúbicos", fontweight='bold')
ax.set_xlabel("Tempo adimensional ($t$)")
ax.set_ylabel("Potência Instantânea $\\mathcal{P}(t)$")
ax.legend(fontsize=11)
ax.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
caminho_salvar_interp = os.path.join(diretorio_atual, "interpolacao_potencia_ex1_ex2.png")
plt.savefig(caminho_salvar_interp, dpi=300)
print(f"✓ Gráfico de interpolação salvo em: {caminho_salvar_interp}")

t_zoom_ini = t_dados[0] + 0.25 * (t_dados[-1] - t_dados[0])
t_zoom_fim = t_dados[0] + 0.50 * (t_dados[-1] - t_dados[0])
mascara_zoom = (t_fino >= t_zoom_ini) & (t_fino <= t_zoom_fim)
P_zoom = np.concatenate([P_cubico[mascara_zoom], P_linear[mascara_zoom]])
margem = 0.1 * (P_zoom.max() - P_zoom.min() + 1e-12)

ax.set_xlim(t_zoom_ini, t_zoom_fim)
ax.set_ylim(P_zoom.min() - margem, P_zoom.max() + margem)
ax.set_title(f"Zoom ({t_zoom_ini:.2f} a {t_zoom_fim:.2f}s) - Analisando Suavidade nas Interfaces",
             fontweight='bold')
caminho_salvar_zoom = os.path.join(diretorio_atual, "interpolacao_potencia_zoom.png")
plt.savefig(caminho_salvar_zoom, dpi=300)
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
print(" EXERCÍCIO 3 (Seção 6.4.3): Regressão Polinomial Global")
print("-" * 85)

# parte da nat vai aqui

print("\n" + "=" * 80)
print("EXERCÍCIO 3 (Seção 6.4.3) CONCLUÍDO.")
print("=" * 80)


# ========================================================================
# SEÇÃO 6.4.3 - EXERCÍCIO 4
# Análise de Sensibilidade a Ruídos Estocásticos
# ========================================================================

print("\n" + "-"*85)
print(" EXERCÍCIO 4 (Seção 6.4.3): Análise de Sensibilidade a Ruídos Estocásticos")
print("-" * 85)

# parte da nat vai aqui

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

# parte do zuffo vai aqui

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
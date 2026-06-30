import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import time

diretorio_atual = os.path.dirname(os.path.abspath(__file__))
diretorio_pai = os.path.dirname(diretorio_atual)
if diretorio_pai not in sys.path:
    sys.path.insert(0, diretorio_pai)

from functionsGD import Setup_Base_GD, GD_Gera_Condutancias_Nominais, GD_Avalia_Vazao_Falhas, Interpola_Potencia_Linear, Interpola_Potencia_Cubica

print("="*85)
print(" GÊMEO DIGITAL - CAPÍTULO 6")
print("="*85)

# ==============================================================================
# SEÇÃO 6.2: ESTADO NOMINAL DO SISTEMA
# ==============================================================================
print("\n[1/2] Configurando a base e resolvendo o campo térmico (Executado apenas 1x)...")
t_setup = time.time()
base_gd = Setup_Base_GD()
C_T_nom = GD_Gera_Condutancias_Nominais(base_gd)
print(f"✓ Base inicializada em {time.time() - t_setup:.2f} segundos.")


# ==============================================================================
# SEÇÃO 6.3.2 (EXERCÍCIO 1): COMPORTAMENTO ASSINTÓTICO (Prob vs p_fail)
# ==============================================================================
print("\n" + "-"*85)
print(" EXERCÍCIO 1: PROBABILIDADE DE FALHA DA VAZÃO vs P_FAIL")
print("-" * 85)

N_simulacoes = 1000  # Tamanho amostral convergido
q_critico = 1.25e-5 
p_inlet = 5000.0    

# Domínio exigido pelo PDF: p_fail em [0.05, 0.65]
p_fails = np.linspace(0.05, 0.65, 13)

# Dois fatores de severidade de obstrução distintos exigidos
fatores_obstrucao = [5.0, 10.0]

# Dicionário para guardar as duas curvas do gráfico
resultados_prob = {5.0: [], 10.0: []}

print(f"Iniciando varredura estocástica (Isso executará {len(p_fails) * 2 * N_simulacoes} simulações rápidas)...")
t0 = time.time()

for f_obs in fatores_obstrucao:
    print(f"\n-> Analisando fator de severidade f_obs = {f_obs}")
    for p_fail in p_fails:
        falhas_acumuladas = 0
        
        # Roda o Monte Carlo N vezes para este ponto específico do gráfico
        for _ in range(N_simulacoes):
            q_in = GD_Avalia_Vazao_Falhas(base_gd, C_T_nom, p_inlet, p_fail, f_obs)
            if q_in < q_critico:
                falhas_acumuladas += 1
                
        prob_convergida = falhas_acumuladas / N_simulacoes
        resultados_prob[f_obs].append(prob_convergida)
        print(f"   p_fail = {p_fail:.2f} | Probabilidade = {prob_convergida:.2%}")

print(f"\n✓ Varredura estocástica concluída em {time.time() - t0:.2f} segundos.")

# --- Renderização Gráfica do 3º Bullet Point ---
fig, ax = plt.subplots(figsize=(8, 5.5))

# Plota as duas curvas com as formatações clássicas
ax.plot(p_fails, resultados_prob[5.0], 's-', color='dodgerblue', linewidth=2, label='$f_{obs} = 5$')
ax.plot(p_fails, resultados_prob[10.0], 'o-', color='darkorange', linewidth=2, label='$f_{obs} = 10$')

ax.set_title("Probabilidade de Falha Global vs Obstrução Individual", fontweight='bold')
ax.set_xlabel("Probabilidade de obstrução individual ($p_{fail}$)")
ax.set_ylabel("Prob. Global Convergida P($Q_{inlet} < 1.25\\times 10^{-5}$)")
ax.grid(True, linestyle='--', alpha=0.6)
ax.legend(fontsize=11)

plt.tight_layout()
caminho_salvar = os.path.join(diretorio_atual, "monte_carlo_p_fail_comparativo.png")
plt.savefig(caminho_salvar, dpi=300)
print(f"✓ Gráfico comparativo salvo em: {caminho_salvar}")
plt.show()

# ==============================================================================
# SEÇÃO 6.3.2 (EXERCÍCIO 2): MONTE CARLO DINÂMICO (ESPAÇO RESERVADO)
# ==============================================================================
# TODO: Implementar posteriormente a varredura temporal para avaliar E_total < 7.0






# ==============================================================================
# SEÇÃO 6.4.3: INTERPOLAÇÃO DE DADOS (EXERCÍCIOS 1 E 2)
# ==============================================================================
print("\n" + "-"*85)
print(" EXERCÍCIOS 1 E 2 (6.4.3): INTERPOLAÇÃO LINEAR E CÚBICA DA POTÊNCIA")
print("-" * 85)

# 1. Definindo o vetor de tempo conforme exigido: t em [0, 4] com dt = 0.05
dt = 0.05
t_dados = np.arange(0, 4.0 + dt, dt)

# 2. Dados da Simulação:
# AVISO: Como o loop dinâmico (6.3.2 Ex 2) é um TODO, vamos gerar uma curva 
# de "dados fictícios" (dummy data) imitando a flutuação e estabilização do sistema.
# Quando a dinâmica real estiver pronta, basta substituir P_dados pelos seus dados reais.
P_dados = 10.0 * (1.0 - np.exp(-1.5 * t_dados)) * (1.0 + 0.1 * np.sin(5 * t_dados))

# 3. Gerando os modelos de interpolação
print("-> Construindo Splines Lineares (Ex 1)...")
interpolador_linear = Interpola_Potencia_Linear(t_dados, P_dados)

print("-> Construindo Splines Cúbicos (Ex 2)...")
interpolador_cubico = Interpola_Potencia_Cubica(t_dados, P_dados)

# 4. Avaliando os modelos em uma malha temporal bem fina para ver a "suavidade"
t_fino = np.linspace(0, 4.0, 500)
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

# Zoom adicional automático para facilitar a comparação visual do "ganho de suavidade"
ax.set_xlim(1.0, 2.0)
ax.set_ylim(8.0, 12.0)
ax.set_title("Zoom (1.0 a 2.0s) - Analisando Suavidade nas Interfaces", fontweight='bold')
caminho_salvar_zoom = os.path.join(diretorio_atual, "interpolacao_potencia_zoom.png")
plt.savefig(caminho_salvar_zoom, dpi=300)
print(f"✓ Gráfico com zoom salvo em: {caminho_salvar_zoom}")

plt.show()
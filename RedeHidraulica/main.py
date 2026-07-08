from functions import (Assembly, GeraGrafo, PlotaRede, SolveNetwork, 
                        CalculoCondutancia, createK, createD, calc_vazao, 
                        calc_potencia, AssemblyVectorC)
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")
import time

# CONSTRUÇÃO DO MODELO FÍSICO

Xno, conec = GeraGrafo(levels=3)
mm_to_m = 0.001
Xno = Xno * mm_to_m
C = AssemblyVectorC(conec, Xno)
n_inlet = 0                 # nó mais a esquerda
n_outlet = len(Xno) - 1     # nó mais a direita

print(f"Rede gerada com {Xno.shape[0]} nós e {conec.shape[0]} canos.\n")

# Teste do funcionamento da Rede Hidrúalica

ps = {str(n_outlet): 0} # pressão atmosférica no último nó (outlet)
Qs = {'1': 1.0e-6} # vazão de 1 mL/s no primeiro nó (inlet)

pressure = SolveNetwork(conec, C, ps=ps, Qs=Qs)
matriz_vazao = calc_vazao(conec, C, pressure)

fig, ax = PlotaRede(conec, Xno, pressure, matriz_vazao, factor_units=mm_to_m)
plt.show()

# EXERCÍCIOS PARA INVESTIGAR O COMPORTAMENTO DO SISTEMA (SEÇÃO 1.4.3)

# ==============================================================================
# OBSERVAÇÃO SOBRE A INDEXAÇÃO DOS NÓS NAS CONDIÇÕES DE CONTORNO (ps e Qs)
# ------------------------------------------------------------------------------
# A função SolveNetwork foi originalmente construída com a instrução interna 
# i = int(node) - 1. Isso foi feito para converter numerações naturais 
# (nós de 1 a N) para os índices nativos do Python (que vão de 0 a N-1).
#
# Logo, passar o index do nó do dicionário foi feito da seguinte forma:
# Nó Alvo (n) = Chave Dicionário (n + 1).
# ==============================================================================

# A: Itens 1 e 2

ps_A = {
    '6': 0.0,       # nó 5 -> pressão atmosférica (saída 1)
    '216': 0.0      # nó 215 -> pressão atmosférica (saída 2)
}

Qs_A = {
    '1': 1.0e-7,    # nó 0 -> injeção de vazão 0.1 mL/s
    '176': 1.0e-6    # nó 175 -> injeção de vazão 1 mL/s
}

pressure_A = SolveNetwork(conec, C, ps=ps_A, Qs=Qs_A)
matriz_vazao_A = calc_vazao(conec, C, pressure_A)

fig, ax = PlotaRede(conec, Xno, pressure_A, matriz_vazao_A, factor_units=mm_to_m)
plt.show()

# B: Item 3

ps_B = {str(n_outlet + 1): 0}

# vazão arbitrária de teste (b^(1))
Q_teste = 1.0e-6
Qs_B = {'1': Q_teste}

# pressão gerada por essa vazão de teste (p^(1))
pressure_teste = SolveNetwork(conec, C, ps=ps_B, Qs=Qs_B)
pressao_no_inlet = pressure_teste[0] # pressão no nó de entrada

# escalar alpha para atingir os 100 Pa desejados
pressao_alvo = 100.0
alpha = pressao_alvo / pressao_no_inlet

vazao_inlet_real = alpha * Q_teste

print(f"Escalar alpha = {alpha:.2f}:")
print(f"Vazão necessária para manter 100 Pa no Inlet: {vazao_inlet_real:.4e} m³/s\n")

# C: Itens 4 e 5

t_array = np.linspace(0, 10, 1000)
omega_sin = 3.0
omega_cos = 4.0

ps_base = {str(n_outlet + 1): 0}

# Calcular vetores base unitários (1 m³/s de injeção) para permitir escalonamento
Qs_unit_sin = {'1': 1.0}
p_base_vetor_sin = SolveNetwork(conec, C, ps=ps_base, Qs=Qs_unit_sin)

Qs_unit_cos = {'176': 1.0}
p_base_vetor_cos = SolveNetwork(conec, C, ps=ps_base, Qs=Qs_unit_cos)

pressao_maxima_ex4 = []
pressao_maxima_ex5 = []

for t in t_array:
    # Calcular vazões reais no instante t (em m^3/s)
    # 1 mL/s = 1.0e-6 m^3/s
    Q_t_sin = (1.0 + 0.1 * np.sin(omega_sin * t)) * 1.0e-6
    Q_t_cos = (0.1 + 0.01 * np.cos(omega_cos * t)) * 1.0e-6
    
    # EXERCÍCIO 4: Apenas injeção no nó 0
    p_t_ex4 = p_base_vetor_sin * Q_t_sin
    pressao_maxima_ex4.append(np.max(p_t_ex4))
    
    # EXERCÍCIO 5: Superposição vetorial dos campos de pressão (nós 0 e 175)
    p_t_ex5 = p_base_vetor_sin * Q_t_sin + p_base_vetor_cos * Q_t_cos
    pressao_maxima_ex5.append(np.max(p_t_ex5))

# Plots
plt.figure(figsize=(10, 6))
plt.plot(t_array, pressao_maxima_ex4, color='blue', linestyle='--', linewidth=2, label='Ex 4 (Apenas Nó 0)')
plt.plot(t_array, pressao_maxima_ex5, color='purple', linewidth=2, label='Ex 5 (Nós 0 e 175)')
plt.title("Pressão Máxima na Rede")
plt.xlabel("Tempo (s)")
plt.ylabel("Pressão Máxima na Rede (Pa)")
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc="upper right")
plt.show()

# D: Item 6

t_array = np.linspace(0, 10, 100)

ps_D = {str(n_outlet + 1): 0}
Qs_D = {'1': 1.0e-7}

pressure_base_D = SolveNetwork(conec, C, ps=ps_D, Qs=Qs_D)
max_pressure_base = np.max(pressure_base_D)

max_pressure = []
mu_0 = 0.001  # Viscosidade base a 20°C

for t in t_array:
    T_t = 20 + 0.9 * (t**2)
    mu_t = 0.001791 / (1 + 0.03368 * T_t + 0.000221 * (T_t**2))

    p_max_t = max_pressure_base * (mu_t / mu_0)
    max_pressure.append(p_max_t)

plt.figure(figsize=(10, 5))
plt.plot(t_array, max_pressure, color='red', linewidth=2)
plt.title("Efeito do Aquecimento na Pressão Máxima")
plt.xlabel("Tempo (s)")
plt.ylabel("Pressão Máxima na Rede (Pa)")
plt.grid(True, linestyle='--', alpha=0.7)
plt.show()

# E: Item 7

print(f"{'Nível':<7} | {'Qtd. Nós':<10} | {'T. Montagem (s)':<18} | {'T. Resolução (s)':<18}")
print("-" * 62)

for level in [1,2,3,4]:
    Xno_test, conec_test = GeraGrafo(levels=level)
    Xno_test = Xno_test * mm_to_m
    n_nos = Xno_test.shape[0] #número de nós

    tempos_montagem = []
    tempos_resolucao = []

    for _ in range(10):

        #MONTAGEM:
        t_inicio_montagem = time.perf_counter()

        C_test = AssemblyVectorC(conec_test, Xno_test)
        A_test = Assembly(conec_test, C_test)

        t_fim_montagem = time.perf_counter()
        tempos_montagem.append(t_fim_montagem - t_inicio_montagem)

        #RESOLUÇÃO:

        #Preparação do sistema linear para que não interfira no tempo da resolução
        Atilde_test = A_test.copy()
        idx_out = n_nos - 1 #numero do nó de saída
        Atilde_test[idx_out, :] = 0
        Atilde_test[idx_out, idx_out] = 1
        b_test = np.zeros(n_nos) #Ax = b, onde, no caso, b seria o vetor de vazões
        b_test[0] = 1.0e-7 # injeção de 0.1 mL/s no inlet

        #Medição do tempo
        t_inicio_resolucao = time.perf_counter()
        
        pressure_test = np.linalg.solve(Atilde_test, b_test)

        t_fim_resolucao = time.perf_counter()
        tempos_resolucao.append(t_fim_resolucao - t_inicio_resolucao)

    media_montagem = np.mean(tempos_montagem)
    media_resolucao = np.mean(tempos_resolucao)

    print(f"{level:<7} | {n_nos:<10} | {media_montagem:<18.6f} | {media_resolucao:<18.6f}")

plt.show()
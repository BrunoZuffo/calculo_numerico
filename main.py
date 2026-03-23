import matplotlib
matplotlib.use("TkAgg")

import time

from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, createK, createD, calc_vazao, calc_potencia, AssemblyVectorC
import numpy as np
import matplotlib.pyplot as plt

# CONSTRUÇÃO DO MODELO FÍSICO

Xno, conec = GeraGrafo(levels=3)
mm_to_m = 0.001
Xno = Xno * mm_to_m
C = AssemblyVectorC(conec, Xno)
n_inlet = 0                 # nó mais a esquerda
n_outlet = len(Xno) - 1     # nó mais a direita

print(f"Rede gerada com {Xno.shape[0]} nós e {conec.shape[0]} canos.\n")

# EXERCÍCIOS PARA INVESTIGAR O COMPORTAMENTO DO SISTEMA (SEÇÃO 1.4.3)

# ==============================================================================
# OBSERVAÇÃO SOBRE A INDEXAÇÃO DOS NÓS NAS CONDIÇÕES DE CONTORNO (ps e Qs)
# ------------------------------------------------------------------------------
# A função `SolveNetwork` foi originalmente construída com a instrução interna 
# `i = int(node) - 1`. Isso foi feito para converter numerações amigáveis 
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

print(f"Escalar alpha = {alpha:.2f}):")
print(f"Vazão necessária para manter 100 Pa no Inlet: {vazao_inlet_real:.4e} m³/s\n")

# C: Itens 4 e 5

t_array = np.linspace(0,10,1000) #discretização do tempo: 1000 passos entre 0 e 10s
omega_sin = 3.0 #frequência angular
omega_cos = 4.0

ps_base = {str(n_outlet + 1): 0}

Qs_base_sin = {'1': 1.0e-6} #vazão base do nó de entrada 0 (1 - 1 = 0), 1ml/s
pressao_base_sin = SolveNetwork(conec, C, ps=ps_base, Qs=Qs_base_sin)
pressao_maxima_base_sin = np.max(pressao_base_sin)

Qs_base_cos = {'176': 1.0e-6} #vazão base do nó 175 (176 - 1 = 175), 1ml/s
pressao_base_cos = SolveNetwork(conec, C, ps=ps_base, Qs=Qs_base_cos)
pressao_maxima_base_cos = np.max(pressao_base_cos)

pressao_maxima_tempo_sin = []
pressao_maxima_tempo_cos = []

for t in t_array:
    fator_tempo_sin = 1.0 + 0.1 * np.sin(omega_sin * t)
    fator_tempo_cos = 0.1 + 0.01 * np.cos(omega_cos * t)
    
    p_max_atual_sin = pressao_maxima_base_sin * fator_tempo_sin
    p_max_atual_cos = pressao_maxima_base_cos * fator_tempo_cos
    pressao_maxima_tempo_sin.append(p_max_atual_sin)
    pressao_maxima_tempo_cos.append(p_max_atual_cos)

plt.figure(figsize=(10, 6))
plt.plot(t_array, pressao_maxima_tempo_sin, color='blue', linewidth=2, label=f'Injeção com sin({omega_sin}t)')
plt.plot(t_array, pressao_maxima_tempo_cos, color='orange', linewidth=2, label=f'Injeção com cos({omega_cos}t)')

plt.title("Pressão Máxima com Injeção Dinâmica (Itens 4 e 5)")
plt.xlabel("Tempo (s)")
plt.ylabel("Pressão Máxima na Rede (Pa)")
plt.grid(True, linestyle=':', alpha=0.7)
plt.legend(loc="upper right")

# D: Item 6

t_array = np.linspace(0,10,100) #discretização do tempo: 100 passos entre 0 e 10s

ps_D = {
    str(n_outlet + 1): 0
}
Qs_D = {
    '1': 1.0e-7 #0.1 mL/s
}

max_pressure = [] #vetor para guardar os valores máximos de pressão ao longo do tempo

for t in t_array:
    T_t = 20 + 0.9*(t**2)
    mu_t = 0.001791 / (1 + 0.03368 * T_t + 0.000221 * (T_t**2))

    #Calculo da condutância:
        #a condutância C já calculada depende da geometria do cano e da viscosidade da água a 20 graus celcius (0.001). Para calcular a condutância em função de 
        #mu_t basta multiplicar C por 0.001 (para eliminar o uso da viscosidade a 20 graus) e dividir por mu_t (para passar a usar a viscosidade em função do tempo)
    C_t = C * (0.001/mu_t)

    pressure_t = SolveNetwork(conec, C_t, ps=ps_D, Qs=Qs_D)
    max_pressure.append(np.max(pressure_t))

plt.figure(figsize=(10, 5))
plt.plot(t_array, max_pressure, color='red', linewidth=2)
plt.title("Efeito do Aquecimento na Pressão Máxima (Item 6)")
plt.xlabel("Tempo (s)")
plt.ylabel("Pressão Máxima na Rede (Pa)")
plt.grid(True, linestyle='--', alpha=0.7)

# E: Item 7

#Cabeçalho da tabela
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

    #Linha da tabela
    print(f"{level:<7} | {n_nos:<10} | {media_montagem:<18.6f} | {media_resolucao:<18.6f}")


# RESULTADOS E PLOTS

#Rede Hidráulica
fig, ax = PlotaRede(conec, Xno, pressure_A, matriz_vazao_A, factor_units=mm_to_m)

plt.show()
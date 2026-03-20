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

# A: Itens 1 e 2

ps_A = {
    '1': 100,   # nó 1 → alta pressão (entrada)
    '4': 0      # nó 4 → baixa pressão (saída)
}

Qs_A = {
    '3': 10,     # nó 3 → injeção de vazão 10
    '100': 200   # nó 100 → injeção de vazão 200
}

pressure_A = SolveNetwork(conec, C, ps=ps_A, Qs=Qs_A)
matriz_vazao_A = calc_vazao(conec, C, pressure_A)

# B: Item 3

ps_B = {
    '1': 100,
    str(n_outlet + 1): 0 #"+ 1" por conta da subtração em SolveNetwork
}

Qs_B = {} #vazio pois não há injeção de vazão

pressure_B = SolveNetwork(conec, C, ps=ps_B, Qs=Qs_B)
A_original = Assembly(conec, C) #recupera a matriz original (sem ser substituida pelas linhas com 0 e valores de pressão)
Q_externas = A_original @ pressure_B #calcula o vetor de vazões externas pelo método Ap = Q_ext, agora sabendo o vetor p (pressão recém calculada)
vazao_inlet = Q_externas[0] #vazão do inlet = vazão do nó 0

print(f"Vazão entrando pelo Inlet (para manter pressão 100): {vazao_inlet:.4e} m³/s\n")

# C: Itens 4 e 5

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

#Item 6
plt.figure(figsize=(10, 5))
plt.plot(t_array, max_pressure, color='red', linewidth=2)
plt.title("Efeito do Aquecimento na Pressão Máxima (Item 6)")
plt.xlabel("Tempo (s)")
plt.ylabel("Pressão Máxima na Rede (Pa)")
plt.grid(True, linestyle='--', alpha=0.7)

#Rede Hidráulica
fig, ax = PlotaRede(conec, Xno, pressure_A, matriz_vazao_A, factor_units=mm_to_m)

plt.show()
from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, CalculoCondutancia  
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# matriz de conectividade do exemplo
conec = np.array([
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 2],
    [5, 3],
    [5, 1]
])


C = np.array([2, 2, 1, 2, 1, 2, 2]) # valores de Ck do exemplo

matriz = Assembly(conec, C)

print(matriz)
Qs_teste = {1:3}

matriz = SolveNetwork(conec,C,3,Qs_teste)

print(matriz)



print("1. Gerando a rede hidráulica")
Xno, conec_bruto = GeraGrafo(levels=3)
conec = conec_bruto + 1 

print("2. Calculando as condutâncias dos canos")
C = []
nc = len(conec)
for k in range(nc):
    n1 = conec[k, 0] - 1
    n2 = conec[k, 1] - 1
    
    x1, y1 = Xno[n1]
    x2, y2 = Xno[n2]
    
    # Comprimento do cano
    Lk = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    Ck = CalculoCondutancia(Lk)
    C.append(Ck)

print("3. Preparando as condições de contorno...")
no_entrada = 0
natm = len(Xno) # Outlet no último nó
Qs_base = {no_entrada: 1*1e-6} 

print("4. Resolvendo o sistema base...")
pressao_base = SolveNetwork(conec, C, natm, Qs_base)
pressao_maxima_base = np.max(pressao_base)

print("5. Calculando a simulação no tempo...")
omega = 3.0
tempos = np.linspace(0.0, 10.0, 1000)
pressao_maxima_tempo = []


for t in tempos:
    
    fator_tempo =  1.0 + 0.1 * np.sin(omega * t) 
    p_max_atual = pressao_maxima_base * fator_tempo
    pressao_maxima_tempo.append(p_max_atual)

print("6. Gerando o gráfico da Questão 4")
plt.figure(figsize=(10, 6))
plt.plot(tempos, pressao_maxima_tempo, color='b', label='Pressão Máxima na Rede')
plt.title('Comportamento da Pressão Máxima ao longo do Tempo')
plt.xlabel('Tempo (s)')
plt.ylabel('Pressão máxima (p)')
plt.grid(True)
plt.legend()
plt.show()

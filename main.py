from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, CalculoCondutancia  
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")

from functions import Assembly, GeraGrafo, PlotaRede, SolveNetwork, createK, createD, calc_vazao, calc_potencia, AssemblyVectorC
import numpy as np
import matplotlib.pyplot as plt

## matriz de conectividade do exemplo
#conec = np.array([
#    [1, 2],
#    [2, 3],
#    [3, 4],
#    [4, 5],
#    [5, 2],
#    [5, 3],
#    [5, 1]
#])
#
##coordenadas dos nós
#Xno = np.array([
#    [0, 0],   # nó 1
#    [1, 0],   # nó 2
#    [2, 0],   # nó 3
#    [3, 0],   # nó 4
#    [1.5, 1]  # nó 5
#])
#
## valores de Ck do exemplo
#C = np.array([2, 2, 1, 2, 1, 2, 2]) 

#teste para exercicio 2
ps = {
    1: 100,   # nó 1 → alta pressão (entrada)
    4: 0      # nó 4 → baixa pressão (saída)
}

Qs = {
    3: 10,     # nó 3 → injeção de vazão
    100: 200
}

Xno, conec = GeraGrafo(levels=3)

mm_to_m = 0.001
Xno = Xno * mm_to_m

C = AssemblyVectorC(conec, Xno)

natm = len(Xno) - 1  #nó atmosférico (= 3 no exemplo numerico)
nbomba = 0  #nó conectado à bomba
Qbomba = 1.0e-7  #vazao da bomba (= 3 no exemplo numerico)

Qs = {nbomba: Qbomba} # dicionario

print("Nós:", Xno.shape[0])
print("Conexões:", conec.shape[0])

matriz = Assembly(conec, C)

print("MATRIZ A")
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
pressure = SolveNetwork(conec,C, ps=ps, Qs=Qs)

print("PRESSURE:")
print(pressure)

matriz_vazao = calc_vazao(conec, C, pressure)
print("MATRIZ DE VAZÃO:")
print (matriz_vazao)

print("coord shape:", Xno.shape)
print("conec max index:", conec.max())

fig, ax = PlotaRede(conec, Xno, pressure, matriz_vazao, factor_units=mm_to_m)

plt.show()

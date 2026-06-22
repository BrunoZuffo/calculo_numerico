import numpy as np
import matplotlib.pyplot as plt
from functionsGD import ResolverGemeoDigital
from ..RedeHidraulica import GeraGrafo, createD

# =============================================================================
# 1. ESPECIFICAÇÕES NOMINAIS DO GÊMEO DIGITAL (Ref: Seção 6.2)
# =============================================================================

# Constantes geométricas
Lx = 0.03   # 3 cm
Ly = 0.015  # 1.5 cm

# Dicionário da Placa Térmica
params_placa = {
    'Lx': Lx,
    'Ly': Ly,
    'k_cond': 0.25,        # Condutividade Térmica W/(K m)
    'Q_fonte': 5e5,        # Geração de calor W/m^3
    'T_L': 10.0,           # Borda Esquerda (Celsius)
    'T_R': 30.0,           # Borda Direita (Celsius)
    # Condições de contorno dependentes da posição x
    'T_B_func': lambda x: 10 + 20 * (x / Lx),
    'T_T_func': lambda x: 10 + 20 * (x / Lx),
    
    # Região interna restrita circular
    'circ_xc': 0.0225,     # 2.25 cm
    'circ_yc': 0.0075,     # 0.75 cm
    'circ_raio': 0.0025,   # 0.25 cm
    'circ_T': 35.0,        # Fixo em 35 Celsius
    
    # Parâmetros numéricos da malha da placa
    'nx': 150, 
    'ny': 75
}

# Gera o grafo de Nível 3 (como você fazia no seu Cap 3)
Xno, conec = GeraGrafo(levels=3)
mm_to_m = 0.001
Xno = Xno * mm_to_m  # Convertendo para metros

n_nos = len(Xno)
n_inlet = 0
n_outlet = n_nos - 1

# Matriz de incidência (Usando a sua função createD)
# Nota: Se o seu createD exigir apenas o 'conec', apague o ', n_nos'
matriz_incidencia = createD(conec, n_nos) 

# Coordenadas centrais de cada canal
# O Gêmeo Digital precisa saber o "centro" de cada canal para pegar a temperatura da placa
# O centro é simplesmente a média entre as coordenadas (x,y) dos dois nós do canal
n_canais = len(conec)
coords_canais = np.zeros((n_canais, 2))

for i in range(n_canais):
    # Pega os índices dos nós inicial e final do canal 'i'
    n1, n2 = int(conec[i, 0]), int(conec[i, 1])
    
    # Faz a média do X e a média do Y
    coords_canais[i, 0] = (Xno[n1, 0] + Xno[n2, 0]) / 2.0
    coords_canais[i, 1] = (Xno[n1, 1] + Xno[n2, 1]) / 2.0

# Índices dos canais ligados ao Inlet (Nó 0)
# Busca na matriz 'conec' quais linhas possuem o nó 0
indices_canais_inlet = np.where((conec[:, 0] == n_inlet) | (conec[:, 1] == n_inlet))[0].tolist()

# Dicionário da Rede Hidráulica
params_rede = {
    'nivel_complexidade': 3,
    'largura_canal': 1000e-6,      # 1000 micrometros em metros
    'p_inlet': 5000.0,             # Pressão no nó 0 em Pascal
    'no_inlet': n_inlet,
    'no_outlet': n_outlet,
    
    # Aqui entram os parâmetros calculados automaticamente ali em cima!
    'matriz_incidencia': matriz_incidencia,
    'coords_canais': coords_canais,
    'indices_canais_inlet': indices_canais_inlet
}

# Dicionário da Membrana Elástica
params_membrana = {
    'R': 0.0025,           # Raio 0.25 cm
    'h_espessura': 0.1e-3, # 0.1 mm
    'rho': 900.0,          # Densidade kg/m^3
    'tau': 200.0,          # Tensão mecânica N/m
    'beta_hat': 0.1,       # Atrito viscoso adimensional
    
    # Parâmetros numéricos para a integração temporal e espacial
    'dt': 1e-4,            # Passo de tempo
    't_final': 0.1,        # Tempo final de observação
    'nr': 51               # Quantidade de nós na malha radial
}


# =============================================================================
# 2. VALIDAÇÃO NOMINAL (SINGLE RUN)
# =============================================================================

if __name__ == "__main__":
    
    print("Iniciando a validação do Gêmeo Digital (Condição Nominal sem falhas)...")
    
    # Substitua os "None" no params_rede pelas matrizes geradas pela sua malha de Nível 3
    # Exemplo: params_rede['matriz_incidencia'] = np.loadtxt('matriz_nivel3.csv')
    
    # Execução única
    resultados = ResolverGemeoDigital(
        params_placa=params_placa,
        params_rede=params_rede,
        params_membrana=params_membrana,
        falha_kwargs=None # Estado Ideal, sem falhas
    )
    
    # Coletando os QoIs nominais para baseline
    q_in_nom = resultados['q_inlet_total']
    p_out_nom = resultados['p_outlet']
    w_max_nom = resultados['w_max_absoluto']
    
    print("\n--- RESULTADOS NOMINAIS ---")
    print(f"Vazão de Entrada (Q_in): {q_in_nom:.4e} m³/s")
    print(f"Pressão no Outlet (p_out): {p_out_nom:.2f} Pa")
    print(f"Deflexão Máxima da Membrana (W_max): {w_max_nom:.4e} m")
    
    # Plotagem simples para validar a dinâmica da membrana
    plt.figure(figsize=(8, 5))
    plt.plot(resultados['tempo'], resultados['w_centro'], color='blue', linewidth=2)
    plt.title('Dinâmica Nominal da Membrana (Acoplada)')
    plt.xlabel('Tempo (s)')
    plt.ylabel('Deflexão no Centro (m)')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    print("\nO Gêmeo Digital está montado. A arquitetura está pronta para receber o loop de Monte Carlo no main.")
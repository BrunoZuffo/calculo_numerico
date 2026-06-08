import numpy as np

def TemperaturaMicrocanais(Xno, conec, T_grid, Lx, Ly):
    """
    Função para o Exercício 1 (Seção 4.3.3).
    Dada a distribuição de temperatura T_grid e o grafo da rede (Xno, conec),
    identifica a temperatura local no ponto médio de cada microcanal.
    """
    nc = len(conec)            # Número de microcanais
    T_canais = np.zeros(nc)    # Vetor para guardar a temperatura de cada canal
    
    # Pegamos as dimensões da malha térmica diretamente do formato de T_grid
    Ny, Nx = T_grid.shape
    
    # Os passos espaciais da grelha térmica
    hx = Lx / (Nx - 1)
    hy = Ly / (Ny - 1)
    
    for k in range(nc):
        # 1. Pega os nós de início (n1) e fim (n2) do canal k
        n1 = int(conec[k, 0])
        n2 = int(conec[k, 1])
        
        # 2. Pega as coordenadas reais desses nós
        x1, y1 = Xno[n1, 0], Xno[n1, 1]
        x2, y2 = Xno[n2, 0], Xno[n2, 1]
        
        # 3. Calcula as coordenadas do ponto médio do canal
        xm = (x1 + x2) / 2.0
        ym = (y1 + y2) / 2.0
        
        # 4. Converte as coordenadas físicas (xm, ym) em índices da matriz T_grid
        # Usa o round para pegar o nó da malha térmica que está mais perto do ponto médio
        i_idx = int(round(xm / hx))
        j_idx = int(round(ym / hy))
        
        # Garante que os índices não saem dos limites da matriz devido a arredondamentos
        i_idx = max(0, min(i_idx, Nx - 1))
        j_idx = max(0, min(j_idx, Ny - 1))
        
        # 5. Guarda a temperatura encontrada para esse canal
        T_canais[k] = T_grid[j_idx, i_idx]
        
    return T_canais
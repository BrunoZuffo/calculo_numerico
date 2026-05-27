import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import time

# =====================================================================
# FUNÇÕES OBRIGATÓRIAS (FORNECIDAS PELO PROFESSOR)
# =====================================================================

def calcular_distancia_ponto_segmento(p, a, b):
    """Calcula a menor distância entre um ponto 'p' e o segmento 'ab'."""
    ab = b - a
    ap = p - a
    ab_2 = np.dot(ab, ab)
    if ab_2 == 0:
        return np.linalg.norm(p - a)
    t = np.dot(ap, ab) / ab_2
    t = np.clip(t, 0.0, 1.0)
    ponto_proximo = a + t * ab
    return np.linalg.norm(p - ponto_proximo)

def CreateMapDistance(Lx, Ly, Nx, Ny, nos_rede, conexoes_rede, d_max):
    """
    Mapeia as arestas próximas a cada ponto da grade utilizando cKDTree.
    """
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
    pontos_grade = np.column_stack((X.ravel(), Y.ravel()))
    
    arvore_grade = cKDTree(pontos_grade)
    mapa_proximidade = [[] for _ in range(Nx * Ny)]
    
    for id_aresta, (idx_i, idx_j) in enumerate(conexoes_rede):
        p_i = nos_rede[idx_i]
        p_j = nos_rede[idx_j]
        
        ponto_medio = (p_i + p_j) / 2.0
        comprimento_aresta = np.linalg.norm(p_j - p_i)
        raio_busca = (comprimento_aresta / 2.0) + d_max
        
        indices_candidatos = arvore_grade.query_ball_point(ponto_medio, raio_busca)
        
        for idx_ponto in indices_candidatos:
            ponto = pontos_grade[idx_ponto]
            dist = calcular_distancia_ponto_segmento(ponto, p_i, p_j)
            
            if dist <= d_max:
                mapa_proximidade[idx_ponto].append((id_aresta, dist))
                
    return mapa_proximidade

# =====================================================================
# ROTINAS FÍSICO-MATEMÁTICAS
# =====================================================================

def ObterCondutividadeFaces_ViaNos(Nx, Ny, Lx, Ly, nos_rede, conexoes_rede, d_max, k0):
    # 1. Utiliza a função obrigatória para obter as distâncias nos nós (centro dos volumes)
    mapa_proximidade = CreateMapDistance(Lx, Ly, Nx, Ny, nos_rede, conexoes_rede, d_max)
    
    # 2. Calcula k nos nós da malha 2D
    K_nos_2D = np.zeros((Nx, Ny))
    for j in range(Ny):
        for i in range(Nx):
            idx = i + j * Nx 
            vizinhos = mapa_proximidade[idx]
            if vizinhos:
                soma = sum(1.0 / (1.0 + d) for _, d in vizinhos)
                K_nos_2D[i, j] = k0 * (1.0 + soma)
            else:
                K_nos_2D[i, j] = k0
                
    # 3. Interpola os valores nodais para as faces leste e norte (média aritmética)
    k_e = np.zeros((Nx - 1, Ny))
    for i in range(Nx - 1):
        for j in range(Ny):
            k_e[i, j] = (K_nos_2D[i, j] + K_nos_2D[i+1, j]) / 2.0
            
    k_n = np.zeros((Nx, Ny - 1))
    for i in range(Nx):
        for j in range(Ny - 1):
            k_n[i, j] = (K_nos_2D[i, j] + K_nos_2D[i, j+1]) / 2.0
            
    return k_e, k_n

def CriarSistemaSolidoCondutividadeVariavel(Nx, Ny, Lx, Ly, k_e, k_n, TL, TR, TB, TT, fonte_calor, R_incl, xincl, yincl, TC):
    nunk = Nx * Ny
    A = np.zeros((nunk, nunk))
    b = np.zeros(nunk)
    hx = Lx / (Nx - 1)
    hy = Ly / (Ny - 1)
    
    x_coords = np.linspace(0, Lx, Nx)
    y_coords = np.linspace(0, Ly, Ny)
    
    for i in range(Nx):
        for j in range(Ny):
            Ic = i + j * Nx
            xc, yc = x_coords[i], y_coords[j]
            
            if (xc - xincl)**2 + (yc - yincl)**2 <= R_incl**2:
                A[Ic, Ic] = 1.0
                b[Ic] = TC
                continue
                
            if i == 0:
                A[Ic, Ic] = 1.0; b[Ic] = TL
            elif i == Nx - 1:
                A[Ic, Ic] = 1.0; b[Ic] = TR
            elif j == 0:
                A[Ic, Ic] = 1.0; b[Ic] = TB[i]
            elif j == Ny - 1:
                A[Ic, Ic] = 1.0; b[Ic] = TT[i]
            else:
                Ie = (i + 1) + j * Nx
                Iw = (i - 1) + j * Nx
                In = i + (j + 1) * Nx
                Is = i + (j - 1) * Nx
                
                ke = k_e[i, j]      
                kw = k_e[i-1, j]    
                kn = k_n[i, j]      
                ks = k_n[i, j-1]    
                
                A[Ic, Ic] = (ke + kw) / hx**2 + (kn + ks) / hy**2
                A[Ic, Ie] = -ke / hx**2
                A[Ic, Iw] = -kw / hx**2
                A[Ic, In] = -kn / hy**2
                A[Ic, Is] = -ks / hy**2
                b[Ic] = fonte_calor
                
    return A, b

def calcular_temperatura_media_arestas(conec, Xno, T_solid_flat, Nx, Ny, Lx, Ly, num_subintervalos=100):
    t_inicio = time.perf_counter()
    
    x_coords = np.linspace(0, Lx, Nx)
    y_coords = np.linspace(0, Ly, Ny)
    T_2D = T_solid_flat.reshape((Ny, Nx)).T 
    
    interp = RegularGridInterpolator((x_coords, y_coords), T_2D, method='linear')
    
    nc = len(conec)
    T_arestas = np.zeros(nc)
    
    t = np.linspace(0, 1, num_subintervalos + 1)
    dt = 1.0 / num_subintervalos
    
    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        p1, p2 = Xno[n1], Xno[n2]
        
        pts = p1 + t[:, np.newaxis] * (p2 - p1)
        T_amostras = interp(pts)
        
        integral = (T_amostras[0] + T_amostras[-1]) / 2.0
        if num_subintervalos > 1:
            integral += np.sum(T_amostras[1:-1])
            
        T_arestas[k] = integral * dt
        
    t_fim = time.perf_counter()
    return T_arestas, (t_fim - t_inicio)

def calcular_temperatura_media_arestas_ponto_medio(conec, Xno, T_solid_flat, Nx, Ny, Lx, Ly, num_subintervalos=100):
    t_inicio = time.perf_counter()
    x_coords = np.linspace(0, Lx, Nx)
    y_coords = np.linspace(0, Ly, Ny)
    T_2D = T_solid_flat.reshape((Ny, Nx)).T 
    
    interp = RegularGridInterpolator((x_coords, y_coords), T_2D, method='linear')
    nc = len(conec)
    T_arestas = np.zeros(nc)
    
    dt = 1.0 / num_subintervalos
    # Ponto médio avalia exatamente no meio de cada subintervalo
    s = np.linspace(dt/2.0, 1.0 - dt/2.0, num_subintervalos)
    
    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        p1, p2 = Xno[n1], Xno[n2]
        pts = p1 + s[:, np.newaxis] * (p2 - p1)
        T_amostras = interp(pts)
        T_arestas[k] = np.mean(T_amostras) # A média amostral é equivalente à integral neste domínio
        
    t_fim = time.perf_counter()
    return T_arestas, (t_fim - t_inicio)

def calcular_viscosidade(T):
    return 0.001791 / (1.0 + 0.03368 * T + 0.000221 * T**2)

def atualiza_condutancias(conec, Xno, T_arestas, Area_canal):
    nc = len(conec)
    C_k = np.zeros(nc)
    D_h = np.sqrt(4 * Area_canal / np.pi)
    
    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        Lk = np.linalg.norm(Xno[n2] - Xno[n1])
        mu = calcular_viscosidade(T_arestas[k])
        C_k[k] = (np.pi * D_h**4) / (128 * mu * Lk)
        
    return C_k

# -----------------------------------------------------------------------------
# FUNÇÃO - Plot dos nós com escala de temperatura (Exercício 2, parte 1)
# -----------------------------------------------------------------------------

def plot_nos_cromaticos_hidraulica(conec, Xno, T_nodes, method, Nx, Ny, ax=None):
    """
    Plota o grafo da rede hidráulica mapeando a temperatura nos nós
    dentro de um eixo específico para evitar quebras de janelas Tkinter.
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    # 1. Plotagem das arestas de suporte em preto fino
    for k in range(len(conec)):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        x_aresta = [Xno[n1, 0], Xno[n2, 0]]
        y_aresta = [Xno[n1, 1], Xno[n2, 1]]
        ax.plot(x_aresta, y_aresta, 'k-', linewidth=0.5, zorder=1, alpha=0.6)

    # 2. Plotagem dos nós coloridos de acordo com a temperatura interpolada
    sc = ax.scatter(
        Xno[:, 0], Xno[:, 1],
        c=T_nodes,
        cmap='jet',
        s=15,
        zorder=2,
        edgecolors='black',
        linewidths=0.2
    )

    ax.set_aspect('equal')
    ax.set_title(f'Nós da Rede - {method.capitalize()}', fontsize=10, fontweight='bold')
    ax.set_xlabel('X (m)', fontsize=8)
    ax.set_ylabel('Y (m)', fontsize=8)
    ax.grid(True, linestyle='--', alpha=0.3)
    
    return sc  # Retorna o objeto scatter para podermos criar a colorbar centralizada

# -----------------------------------------------------------------------------
# FUNÇÃO - Plot dos arestas com escala de temperatura (Exercício 3, parte 1)
# -----------------------------------------------------------------------------

def plot_arestas_cromaticas_hidraulics(conec, Xno, T_arestas, titulo="Temperatura Média nas Arestas"):
    fig, ax = plt.subplots(figsize=(10, 5))

    segmentos = []
    for k in range(len(conec)):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        segmentos.append((Xno[n1], Xno[n2]))

    norm = plt.Normalize(vmin=T_arestas.min(), vmax=T_arestas.max())
    lc = LineCollection(segmentos, cmap='jet', norm=norm, linewidths=2.5, zorder=1)
    lc.set_array(T_arestas)
    ax.add_collection(lc)

    ax.scatter(Xno[:, 0], Xno[:, 1], c='black', s=10, zorder=2)
    
    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label('Temperatura (°C)', fontsize=11)
    
    ax.set_title(titulo, fontsize=12, fontweight='bold')
    ax.set_xlabel('Posição X (m)', fontsize=10)
    ax.set_ylabel('Posição Y (m)', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()

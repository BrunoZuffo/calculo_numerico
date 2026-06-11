import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
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

# =====================================================================
# EXERCÍCIO 4.3.3 - ITEM 2: TERMO FONTE GAUSSIANO DA REDE
# =====================================================================
# =============================================================================
# FUNÇÕES DO EXERCÍCIO 2 (Secção 4.3.2) - ADICIONADAS
# =============================================================================

def normalizar_rede_hidraulica(Xno_local, conec_local, Lx, Ly):
    Xno_local = np.asarray(Xno_local, dtype=float)
    conec_local = np.asarray(conec_local, dtype=int)

    if Xno_local.ndim != 2:
        raise RuntimeError("Xno deveria ser uma matriz bidimensional.")

    if Xno_local.shape[0] == 2 and Xno_local.shape[1] != 2:
        Xno_local = Xno_local.T

    if Xno_local.shape[1] != 2:
        raise RuntimeError("Xno deveria ter duas colunas: x e y.")

    if conec_local.ndim != 2:
        raise RuntimeError("conec deveria ser uma matriz bidimensional.")

    if conec_local.shape[0] == 2 and conec_local.shape[1] != 2:
        conec_local = conec_local.T

    if conec_local.shape[1] != 2:
        raise RuntimeError("conec deveria ter duas colunas.")

    if np.min(conec_local) == 1:
        conec_local = conec_local - 1

    if np.max(np.abs(Xno_local)) > 0.5:
        Xno_local = Xno_local * 0.001

    if np.min(Xno_local[:, 1]) < -1e-12:
        Xno_local[:, 1] = Xno_local[:, 1] + 0.5 * Ly

    if np.max(conec_local) >= Xno_local.shape[0]:
        raise RuntimeError(
            "conec possui índice maior que o número de nós em Xno.\n"
            "Verifique se a conectividade começa em 0 ou 1."
        )

    return Xno_local, conec_local


def carregar_rede_hidraulica(levels=3, spine_length=6, Lx=0.03, Ly=0.015, globals_dict=None):
    if globals_dict is None:
        globals_dict = globals()

    nomes_possiveis = [
        "generate_graph_arrays",
        "gera_grafo_new",
        "GeraGrafo_new",
        "gera_grapho_new",
        "GeraGrapho_new",
        "GeraGrafo",
        "gera_grafo"
    ]

    funcao_grafo = None
    for nome in nomes_possiveis:
        if nome in globals_dict:
            funcao_grafo = globals_dict[nome]
            print(f"Função de grafo encontrada: {nome}")
            break

    if funcao_grafo is not None:
        tentativas = [
            lambda: funcao_grafo(levels=levels),
            lambda: funcao_grafo(levels),
            lambda: funcao_grafo(levels, spine_length),
            lambda: funcao_grafo(complex_level=levels, spine_length=spine_length),
            lambda: funcao_grafo()
        ]

        saida = None
        ultimo_erro = None

        for tentativa in tentativas:
            try:
                saida = tentativa()
                break
            except TypeError as erro:
                ultimo_erro = erro

        if saida is None:
            raise RuntimeError(
                "Não foi possível chamar a função de geração de grafo.\n"
                f"Último erro: {ultimo_erro}"
            )

        if not isinstance(saida, tuple) or len(saida) < 2:
            raise RuntimeError(
                "A função de grafo deveria retornar algo do tipo:\n"
                "Xno, conec"
            )

        A = np.asarray(saida[0])
        B = np.asarray(saida[1])

        if np.issubdtype(A.dtype, np.integer) and not np.issubdtype(B.dtype, np.integer):
            conec_local = A
            Xno_local = B
        else:
            Xno_local = A
            conec_local = B

        return normalizar_rede_hidraulica(Xno_local, conec_local, Lx, Ly)

    if "Xno" in globals_dict and "conec" in globals_dict:
        print("Usando Xno e conec já existentes na memória.")
        return normalizar_rede_hidraulica(globals_dict["Xno"], globals_dict["conec"], Lx, Ly)

    raise RuntimeError(
        "Nenhuma função de geração de grafo foi encontrada.\n"
        "Execute primeiro a célula que define gera_grafo_new, "
        "gera_grapho_new ou generate_graph_arrays."
    )


def dist_pontos_segmento_vetorizado(P, A, B):
    AB = B - A
    AP = P - A
    AB_len2 = np.dot(AB, AB)

    if AB_len2 < 1e-20:
        return np.linalg.norm(P - A, axis=1)

    t = np.sum(AP * AB, axis=1) / AB_len2
    t = np.clip(t, 0.0, 1.0)
    proj = A + t[:, np.newaxis] * AB

    return np.linalg.norm(P - proj, axis=1)


def calc_campo_fonte(Nx, Ny, Lx, Ly, Xno, conec, S0, d_max, I_j_array):
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X_grid, Y_grid = np.meshgrid(x, y)

    pts = np.column_stack((X_grid.ravel(), Y_grid.ravel()))
    arvore_pts = cKDTree(pts)

    A_edges = Xno[conec[:, 0]]
    B_edges = Xno[conec[:, 1]]

    sigma = d_max / 2.0
    nc = conec.shape[0]
    Sp = np.zeros(len(pts))

    for j in range(nc):
        A = A_edges[j]
        B = B_edges[j]

        ponto_medio = 0.5 * (A + B)
        comprimento = np.linalg.norm(B - A)
        raio_busca = 0.5 * comprimento + d_max

        idx_cand = arvore_pts.query_ball_point(ponto_medio, raio_busca)
        if len(idx_cand) == 0:
            continue

        idx_cand = np.asarray(idx_cand, dtype=int)
        P = pts[idx_cand]

        dist = dist_pontos_segmento_vetorizado(P, A, B)
        mask = dist <= d_max

        if np.any(mask):
            idx_ok = idx_cand[mask]
            Sp[idx_ok] += (
                I_j_array[j]
                * np.exp(-(dist[mask] ** 2) / (2.0 * sigma ** 2))
            )

    fonte_2d = S0 * Sp
    print(f"S0 = {S0:.1e} | fonte min = {np.min(fonte_2d):.3e} | fonte max = {np.max(fonte_2d):.3e}")
    return fonte_2d.reshape((Ny, Nx))


def SolveSystem_FonteEspacial(Nx, Ny, Lx, Ly, k0, fonte_nominal, cilindro, TR, fonte_2d):
    nunk = Nx * Ny
    hx = Lx / (Nx - 1)
    hy = Ly / (Ny - 1)

    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)

    ae = k0 * hy / hx
    aw = k0 * hy / hx
    an = k0 * hx / hy
    ass = k0 * hx / hy

    rows = []
    cols = []
    data = []

    b = hx * hy * (fonte_nominal + fonte_2d.ravel())

    R = cilindro['R']
    cx = cilindro['cx']
    cy = cilindro['cy']
    Tc = cilindro['Tc']

    def ij2n_local(i, j, Nx_val):
        return i + j * Nx_val

    for j in range(Ny):
        for i in range(Nx):
            idx = ij2n_local(i, j, Nx)
            xi = x[i]
            yj = y[j]

            eh_borda = (i == 0 or i == Nx - 1 or j == 0 or j == Ny - 1)
            eh_cilindro = ((xi - cx) ** 2 + (yj - cy) ** 2 <= R ** 2)

            if eh_borda:
                rows.append(idx)
                cols.append(idx)
                data.append(1.0)

                if j == 0 or j == Ny - 1:
                    b[idx] = 10.0 + 20.0 * (xi / Lx)
                elif i == 0:
                    b[idx] = 10.0
                elif i == Nx - 1:
                    b[idx] = TR

            elif eh_cilindro:
                rows.append(idx)
                cols.append(idx)
                data.append(1.0)
                b[idx] = Tc
            else:
                rows.append(idx)
                cols.append(idx)
                data.append(ae + aw + an + ass)

                rows.append(idx)
                cols.append(ij2n_local(i + 1, j, Nx))
                data.append(-ae)

                rows.append(idx)
                cols.append(ij2n_local(i - 1, j, Nx))
                data.append(-aw)

                rows.append(idx)
                cols.append(ij2n_local(i, j + 1, Nx))
                data.append(-an)

                rows.append(idx)
                cols.append(ij2n_local(i, j - 1, Nx))
                data.append(-ass)

    A = sparse.coo_matrix((data, (rows, cols)), shape=(nunk, nunk)).tocsr()
    T_1D = spsolve(A, b)
    return T_1D.reshape((Ny, Nx))


def plotar_rede_sobre_ax(ax, Xno, conec):
    for n1, n2 in conec:
        ax.plot(
            [Xno[n1, 0], Xno[n2, 0]],
            [Xno[n1, 1], Xno[n2, 1]],
            color='black', linewidth=0.55, alpha=0.55, zorder=5
        )
    ax.scatter(
        Xno[:, 0], Xno[:, 1],
        s=7, color='black', edgecolors='white', linewidths=0.20, zorder=6
    )


def obter_limites_e_niveis_individuais(Z, n_levels=26):
    vmin = np.nanmin(Z)
    vmax = np.nanmax(Z)
    if np.isclose(vmin, vmax):
        delta = max(1.0, abs(vmin) * 0.05 + 1e-12)
        vmin = vmin - delta
        vmax = vmax + delta
    levels = np.linspace(vmin, vmax, n_levels)
    return vmin, vmax, levels
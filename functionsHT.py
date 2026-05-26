import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
import time
from matplotlib.collections import LineCollection

# =====================================================================
# FUNÇÕES DO PROFESSOR (mapea_pontos_prox.py)
# =====================================================================
def calcular_distancia_ponto_segmento(p, a, b):
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
# FUNÇÕES DE ACOPLAMENTO ATUALIZADAS
# =====================================================================
def CalcularComprimentos(Xno, conec):
    Nc = len(conec)
    L_canos = np.zeros(Nc)
    for k, (n1, n2) in enumerate(conec):
        dx = Xno[n2, 0] - Xno[n1, 0]
        dy = Xno[n2, 1] - Xno[n1, 1]
        L_canos[k] = np.sqrt(dx**2 + dy**2)
    return L_canos

def ObterTemperaturaSolidoNosCanais(conec, T_solid_flat, mapa_proximidade):
    """
    Calcula a temperatura média do sólido ao redor de cada canal
    usando os nós mapeados.
    """
    Nc = len(conec)
    T_s_canais = np.zeros(Nc)
    contagem_nos_por_canal = np.zeros(Nc)
    
    # Soma as temperaturas de todos os nós associados a cada canal
    for id_no, arestas_proximas in enumerate(mapa_proximidade):
        for id_aresta, dist in arestas_proximas:
            T_s_canais[id_aresta] += T_solid_flat[id_no]
            contagem_nos_por_canal[id_aresta] += 1
            
    # Tira a média
    for k in range(Nc):
        if contagem_nos_por_canal[k] > 0:
            T_s_canais[k] /= contagem_nos_por_canal[k]
        else:
            # Segurança: se nenhum nó foi mapeado, pega a T ambiente do momento
            T_s_canais[k] = np.mean(T_solid_flat)
            
    return T_s_canais, contagem_nos_por_canal

def AssemblyTermicoFluido(Xno, conec, Q, L_canos, P_canos, T_s_canais, rho, cp, hc, T_inlet, n_inlet):
    # (Mantida igual ao anterior, pois ela só precisa do T_s_canais que agora é mais preciso)
    Nv = len(Xno)
    Af = np.zeros((Nv, Nv))
    bf = np.zeros(Nv)
    
    for k, (n1, n2) in enumerate(conec):
        qk = Q[k]
        n_in, n_out = (n1, n2) if qk >= 0 else (n2, n1)
        abs_qk = abs(qk)
        Lk = L_canos[k]
        Pk = P_canos[k]
        
        if n_out != n_inlet:
            Af[n_out, n_out] += rho * cp * abs_qk
            Af[n_out, n_in]  += -rho * cp * abs_qk
            Af[n_out, n_in]  += hc * Pk * Lk
            bf[n_out]        += hc * Pk * Lk * T_s_canais[k]
            
    Af[n_inlet, :] = 0.0
    Af[n_inlet, n_inlet] = 1.0
    bf[n_inlet] = T_inlet
    return Af, bf

def ModificarSistemaSolido(A_s_base, b_s_base, conec, Q, Tf, L_canos, P_canos, hc, mapa_proximidade, contagem_nos_por_canal):
    """
    Distribui o termo convectivo do canal por todos os nós da placa
    que estão na sua área de influência.
    """
    A_s = A_s_base.copy().tolil() if sparse.issparse(A_s_base) else A_s_base.copy()
    b_s = b_s_base.copy()
    
    for id_no, arestas_proximas in enumerate(mapa_proximidade):
        for id_aresta, dist in arestas_proximas:
            qk = Q[id_aresta]
            n1, n2 = conec[id_aresta]
            n_in = n1 if qk >= 0 else n2
            
            Tf_canal = Tf[n_in]
            Lk = L_canos[id_aresta]
            Pk = P_canos[id_aresta]
            num_nos_mapeados = contagem_nos_por_canal[id_aresta]
            
            # Conservação de energia: a área de troca do canal é dividida
            # igualmente entre todos os nós do sólido afetados por ele
            if num_nos_mapeados > 0:
                fator_area = (hc * Pk * Lk) / num_nos_mapeados
                
                # Adiciona o termo na diagonal e no vetor de forças implicitamente
                A_s[id_no, id_no] += fator_area
                b_s[id_no]        += fator_area * Tf_canal
                
    return A_s.tocsc() if sparse.issparse(A_s) else A_s

def calcular_temperatura_media_arestas(conec, Xno, T_fluid_nodes, num_subintervalos=1000):
    """
    Calcula a temperatura média em cada aresta (canal) da rede hidráulica
    empregando a Regra do Trapézio (Simples ou Composta).
    """
    t_inicio = time.perf_counter()
    nc = len(conec)
    T_arestas = np.zeros(nc)
    
    # Pontos de amostragem adimensionais ao longo do canal (de 0 a 1)
    t = np.linspace(0, 1, num_subintervalos + 1)
    dt = 1.0 / num_subintervalos
    
    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        
        # Temperatura nas extremidades do canal k (valores conhecidos nos nós)
        T1 = T_fluid_nodes[n1]
        T2 = T_fluid_nodes[n2]
        
        # Interpolação linear da temperatura do fluido ao longo do canal nos pontos de quadratura
        T_amostras = T1 + t * (T2 - T1)
        
        # Aplicação da Regra do Trapézio
        integral = (T_amostras[0] + T_amostras[-1]) / 2.0
        if num_subintervalos > 1:
            integral += np.sum(T_amostras[1:-1])
        integral *= dt
        
        T_arestas[k] = integral
        
    t_fim = time.perf_counter()
    tempo_execucao = t_fim - t_inicio
    
    return T_arestas, tempo_execucao


def plot_arestas_cromaticas_hidraulics(conec, Xno, T_nodes, T_arestas, titulo="Rede Hidráulica"):
    """
    Plota o grafo da rede hidráulica mapeando os valores obtidos de temperatura média 
    em uma escala cromática aplicada diretamente sobre as arestas.
    """
    fig, ax = plt.subplots(figsize=(10, 5))

    # Construção dos segmentos de reta para o LineCollection
    segmentos = []
    for k in range(len(conec)):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])
        segmentos.append((Xno[n1], Xno[n2]))

    # Normalização de cores baseada nos limites térmicos dos nós do fluido
    norm = plt.Normalize(vmin=T_nodes.min(), vmax=T_nodes.max())
    
    # Plotagem das arestas coloridas pela temperatura média obtida na quadratura
    lc = LineCollection(segmentos, cmap='jet', norm=norm, linewidths=2.5, zorder=1)
    lc.set_array(T_arestas)
    ax.add_collection(lc)

    # Plotagem dos nós coloridos pela temperatura pontual do fluido
    sc = ax.scatter(
        Xno[:, 0], Xno[:, 1],
        c=T_nodes,
        cmap='jet',
        norm=norm,
        s=30,
        zorder=2,
        edgecolors='black',
        linewidths=0.5
    )
    
    # Barra de cores lateral e configurações dos eixos
    cbar = fig.colorbar(lc, ax=ax)
    cbar.set_label('Temperatura (°C)', fontsize=11)
    
    ax.set_title(titulo, fontsize=12, fontweight='bold')
    ax.set_xlabel('Posição X (m)', fontsize=10)
    ax.set_ylabel('Posição Y (m)', fontsize=10)
    ax.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.show()



def CalcularK_Face(pts_face, nos_rede, conexoes_rede, d_max, k0):
    """
    Identifica as arestas da rede próximas a cada ponto médio de interface
    e calcula a condutividade perturbada k_f.
    """
    arvore = cKDTree(pts_face)
    k_face = np.ones(len(pts_face)) * k0
    mapa_prox = [[] for _ in range(len(pts_face))]
    
    for id_aresta, (idx_i, idx_j) in enumerate(conexoes_rede):
        p_i = nos_rede[idx_i]
        p_j = nos_rede[idx_j]
        ponto_medio = (p_i + p_j) / 2.0
        comp = np.linalg.norm(p_j - p_i)
        raio = (comp / 2.0) + d_max
        
        # Otimiza a busca encontrando apenas pontos potencialmente próximos
        candidatos = arvore.query_ball_point(ponto_medio, raio)
        for idx_pt in candidatos:
            # Calcula a distância exata ao segmento do canal
            dist = calcular_distancia_ponto_segmento(pts_face[idx_pt], p_i, p_j)
            if dist <= d_max:
                mapa_prox[idx_pt].append(dist)
                
    for i, dists in enumerate(mapa_prox):
        if dists:
            # Fórmula: k_f = k_0 * (1 + sum(1 / (1 + d_j)))
            soma = sum(1.0 / (1.0 + d) for d in dists)
            k_face[i] = k0 * (1.0 + soma)
            
    return k_face

def ObterCondutividadeFaces(Nx, Ny, Lx, Ly, nos_rede, conexoes_rede, d_max, k0):
    """
    Gera as matrizes 2D de condutividade para as faces leste (k_e) e norte (k_n).
    """
    hx = Lx / (Nx - 1)
    hy = Ly / (Ny - 1)

    # Faces Leste: pontos médios entre i e i+1
    x_e = np.linspace(hx/2, Lx - hx/2, Nx - 1)
    y_e = np.linspace(0, Ly, Ny)
    X_e, Y_e = np.meshgrid(x_e, y_e, indexing='ij')
    pts_e = np.column_stack((X_e.ravel(), Y_e.ravel()))
    k_e_flat = CalcularK_Face(pts_e, nos_rede, conexoes_rede, d_max, k0)
    k_e = k_e_flat.reshape((Nx - 1, Ny))

    # Faces Norte: pontos médios entre j e j+1
    x_n = np.linspace(0, Lx, Nx)
    y_n = np.linspace(hy/2, Ly - hy/2, Ny - 1)
    X_n, Y_n = np.meshgrid(x_n, y_n, indexing='ij')
    pts_n = np.column_stack((X_n.ravel(), Y_n.ravel()))
    k_n_flat = CalcularK_Face(pts_n, nos_rede, conexoes_rede, d_max, k0)
    k_n = k_n_flat.reshape((Nx, Ny - 1))

    return k_e, k_n

def CriarSistemaSolidoCondutividadeVariavel(Nx, Ny, Lx, Ly, k_e, k_n, TL, TR, TB, TT, fonte_calor, R_incl, xincl, yincl, TC):
    """
    Monta a matriz do sistema de condução térmica usando k variável nas faces.
    Inclui uma máscara circular onde a temperatura é fixada em TC.
    """
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
            
            # Máscara para a inclusão circular
            if (xc - xincl)**2 + (yc - yincl)**2 <= R_incl**2:
                A[Ic, Ic] = 1.0
                b[Ic] = TC
                continue
                
            # Condições de contorno de Dirichlet nas bordas da placa
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
                
                # Valores de k nas interfaces da célula (i,j)
                ke = k_e[i, j]      
                kw = k_e[i-1, j]    
                kn = k_n[i, j]      
                ks = k_n[i, j-1]    
                
                # Conservação de calor - Diferenças Finitas 2D
                A[Ic, Ic] = (ke + kw) / hx**2 + (kn + ks) / hy**2
                A[Ic, Ie] = -ke / hx**2
                A[Ic, Iw] = -kw / hx**2
                A[Ic, In] = -kn / hy**2
                A[Ic, Is] = -ks / hy**2
                b[Ic] = fonte_calor
                
    return A, b

# =====================================================================
# EXERCÍCIO 4 - CAPÍTULO 4
# Temperatura média nas arestas -> viscosidade -> condutância hidráulica
# =====================================================================

from scipy.interpolate import RegularGridInterpolator


def viscosidade_agua_T(T):
    """
    Viscosidade dinâmica da água em função da temperatura T em °C.
    Retorna mu em Pa.s.
    """
    T = np.asarray(T)
    return 0.001791 / (1.0 + 0.03368*T + 0.000221*T**2)


def criar_interpolador_temperatura(T_flat, Nx, Ny, Lx, Ly, metodo="linear"):
    """
    Cria interpolador bidimensional para o campo de temperaturas da placa.

    No projeto, o índice global é:
        I = i + j*Nx

    Por isso:
        T_flat.reshape((Ny, Nx)).T
    gera uma matriz no formato (Nx, Ny), compatível com (x, y).
    """
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)

    T_matriz = np.array(T_flat).reshape((Ny, Nx)).T

    interpolador = RegularGridInterpolator(
        (x, y),
        T_matriz,
        method=metodo,
        bounds_error=False,
        fill_value=None
    )

    return interpolador, T_matriz


def calcular_temperatura_media_arestas_campo(
    conec,
    Xno,
    interpolador_T,
    metodo="trapezio",
    num_subintervalos=10
):
    """
    Calcula a temperatura média de cada aresta da rede hidráulica
    usando o campo térmico da placa.

    metodo:
        "trapezio"
        "ponto_medio"
    """
    t_inicio = time.perf_counter()

    nc = len(conec)
    T_arestas = np.zeros(nc)

    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])

        p1 = Xno[n1]
        p2 = Xno[n2]

        if metodo == "ponto_medio":
            valores = []

            for m in range(num_subintervalos):
                alpha = (m + 0.5) / num_subintervalos
                ponto = p1 + alpha * (p2 - p1)
                valores.append(interpolador_T([ponto])[0])

            T_arestas[k] = np.mean(valores)

        elif metodo == "trapezio":
            alphas = np.linspace(0.0, 1.0, num_subintervalos + 1)
            valores = []

            for alpha in alphas:
                ponto = p1 + alpha * (p2 - p1)
                valores.append(interpolador_T([ponto])[0])

            valores = np.array(valores)

            media_integral = valores[0]/2.0 + np.sum(valores[1:-1]) + valores[-1]/2.0
            media_integral = media_integral / num_subintervalos

            T_arestas[k] = media_integral

        else:
            raise ValueError("metodo deve ser 'trapezio' ou 'ponto_medio'.")

    tempo = time.perf_counter() - t_inicio

    return T_arestas, tempo


def calcular_condutancias_por_temperatura(conec, Xno, T_arestas, area_canal=2.5e-7):
    """
    Calcula a condutância hidráulica de cada canal usando a temperatura média
    da própria aresta.
    """
    nc = len(conec)
    C = np.zeros(nc)

    D = np.sqrt(4.0 * area_canal / np.pi)

    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])

        p1 = Xno[n1]
        p2 = Xno[n2]

        Lk = np.linalg.norm(p2 - p1)
        mu_k = viscosidade_agua_T(T_arestas[k])

        kappa_k = np.pi * D**4 / (128.0 * mu_k)
        C[k] = kappa_k / Lk

    return C

# =====================================================================
# EXERC?CIO 4 - CAP?TULO 4
# Temperatura m?dia nas arestas -> viscosidade -> condut?ncia hidr?ulica
# =====================================================================

from scipy.interpolate import RegularGridInterpolator


def viscosidade_agua_T(T):
    """
    Viscosidade din?mica da ?gua em fun??o da temperatura T em ?C.
    Retorna mu em Pa.s.
    """
    T = np.asarray(T)
    return 0.001791 / (1.0 + 0.03368*T + 0.000221*T**2)


def criar_interpolador_temperatura(T_flat, Nx, Ny, Lx, Ly, metodo="linear"):
    """
    Cria interpolador bidimensional para o campo de temperaturas da placa.
    O ?ndice global usado no projeto ? I = i + j*Nx.
    """
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)

    T_matriz = np.array(T_flat).reshape((Ny, Nx)).T

    interpolador = RegularGridInterpolator(
        (x, y),
        T_matriz,
        method=metodo,
        bounds_error=False,
        fill_value=None
    )

    return interpolador, T_matriz


def calcular_temperatura_media_arestas_campo(
    conec,
    Xno,
    interpolador_T,
    metodo="trapezio",
    num_subintervalos=10
):
    """
    Calcula a temperatura m?dia de cada aresta da rede hidr?ulica
    usando o campo t?rmico da placa.

    metodo:
        "trapezio"
        "ponto_medio"
    """
    t_inicio = time.perf_counter()

    nc = len(conec)
    T_arestas = np.zeros(nc)

    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])

        p1 = Xno[n1]
        p2 = Xno[n2]

        if metodo == "ponto_medio":
            valores = []

            for m in range(num_subintervalos):
                alpha = (m + 0.5) / num_subintervalos
                ponto = p1 + alpha * (p2 - p1)
                valores.append(interpolador_T([ponto])[0])

            T_arestas[k] = np.mean(valores)

        elif metodo == "trapezio":
            alphas = np.linspace(0.0, 1.0, num_subintervalos + 1)
            valores = []

            for alpha in alphas:
                ponto = p1 + alpha * (p2 - p1)
                valores.append(interpolador_T([ponto])[0])

            valores = np.array(valores)

            media_integral = valores[0]/2.0 + np.sum(valores[1:-1]) + valores[-1]/2.0
            media_integral = media_integral / num_subintervalos

            T_arestas[k] = media_integral

        else:
            raise ValueError("metodo deve ser 'trapezio' ou 'ponto_medio'.")

    tempo = time.perf_counter() - t_inicio

    return T_arestas, tempo


def calcular_condutancias_por_temperatura(conec, Xno, T_arestas, area_canal=2.5e-7):
    """
    Calcula a condut?ncia hidr?ulica de cada canal usando a temperatura m?dia
    da pr?pria aresta.
    """
    nc = len(conec)
    C = np.zeros(nc)

    D = np.sqrt(4.0 * area_canal / np.pi)

    for k in range(nc):
        n1, n2 = int(conec[k, 0]), int(conec[k, 1])

        p1 = Xno[n1]
        p2 = Xno[n2]

        Lk = np.linalg.norm(p2 - p1)
        mu_k = viscosidade_agua_T(T_arestas[k])

        kappa_k = np.pi * D**4 / (128.0 * mu_k)
        C[k] = kappa_k / Lk

    return C

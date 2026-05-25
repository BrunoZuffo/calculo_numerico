import numpy as np
from scipy import sparse
from scipy.spatial import cKDTree

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

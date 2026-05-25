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
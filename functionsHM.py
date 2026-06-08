import numpy as np
from scipy import sparse
import scipy.sparse.linalg as splinalg

def Monta_Matriz_Global_Acoplada(K_mem, M_mem, A_net, N1, N2, Lx_ad, Ly_ad, R_ad, h_ad, dt, beta, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet):
    """
    Monta a matriz global em blocos para o sistema acoplado (Hidráulico-Mecânico).
    """
    nm = N1 * N2
    np_net = A_net.shape[0]
    
    # 1. Escalas de Adimensionalização (Usando valores físicos)
    w0 = 0.01 * R_fisico
    pref = sigma * w0 / (R_fisico**2)
    vref = (w0 / R_fisico) * np.sqrt(sigma / (rho * e_esp))
    
    # Fator de conversão para a matriz A da rede (Garantindo que seja esparsa)
    A_scaled = A_net.copy() * (pref / (vref * R_fisico**2))
    
    # Impondo condição de contorno no Inlet na matriz A (P_inlet prescrito)
    # A conversão para LIL facilita zerar a linha em matrizes esparsas
    A_scaled_lil = A_scaled.tolil()
    A_scaled_lil[n_inlet, :] = 0.0
    A_scaled_lil[n_inlet, n_inlet] = 1.0
    A_scaled_csr = A_scaled_lil.tocsr()

    # 2. Matrizes da Membrana (Adimensionais)
    D = beta * M_mem

    # 3. Matriz de Acoplamento U
    # Identificando os nós livres (dentro do círculo) para acoplar a vazão
    U = np.zeros((np_net, nm))
    uns = np.ones(nm)
    
    xc, yc = Lx_ad / 2.0, Ly_ad / 2.0
    x = np.linspace(0, Lx_ad, N1)
    y = np.linspace(0, Ly_ad, N2)
    
    for j in range(N2):
        for i in range(N1):
            dist_quadrada = (x[i] - xc)**2 + (y[j] - yc)**2
            if dist_quadrada >= R_ad**2:
                # O ponto está restrito, logo não entra na somatória do volume
                Ic = i + j * N1
                uns[Ic] = 0.0
                
    # Atribui os pesos da área da membrana APENAS na linha correspondente à saída (outlet)
    U[n_outlet, :] = uns[:]
    U_sparse = sparse.csr_matrix(U)
    
    # 4. Montagem dos Blocos da Matriz Global
    Iden = sparse.identity(nm, format='csr')
    idt = 1.0 / dt
    
    h2_adim = h_ad**2
    
    # Matriz global baseada na equação (5.3) da documentação
    blocks = [
        [idt * Iden, -Iden, None],
        [K_mem, idt * M_mem + D, -U_sparse.T],
        [None, U_sparse * h2_adim, A_scaled_csr]
    ]
    
    Aglob = sparse.bmat(blocks, format='csr')
    
    return Aglob, U_sparse, pref, vref

def Resolve_Passo_Tempo(Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, np_net, n_inlet):
    """
    Resolve um único passo de tempo (Euler Implícito Acoplado)
    """
    # Lado direito do sistema
    b_w = w_n / dt
    b_v = (M_mem.dot(v_n)) / dt
    
    b_p = np.zeros(np_net)
    b_p[n_inlet] = p_inlet_adim  # Pressão prescrita na entrada
    
    # Vetor global do lado direito concatenado
    b_glob = np.concatenate([b_w, b_v, b_p])
    
    # Resolver o sistema linear
    sol = splinalg.spsolve(Aglob, b_glob)
    
    # Extrair as soluções individuais dos blocos do vetor
    w_next = sol[0:nm]
    v_next = sol[nm:2*nm]
    p_next = sol[2*nm:]
    
    return w_next, v_next, p_next

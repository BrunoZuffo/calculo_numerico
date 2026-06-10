import numpy as np
from scipy import sparse
import scipy.sparse.linalg as splinalg
from scipy.sparse.linalg import eigsh
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
from functions import GeraGrafo, AssemblyVectorC, Assembly
from functionsM import BuildMatrizes_Eigen_Circular

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

def roda_exercicio_4(Nx, Ny, h_ad, dt, mu, R_fisico, e_esp, sigma, rho, Lx_ad, Ly_ad, R_ad, sigma_ad, rho_ad, e_ad):
    print("\n--- Iniciando Exercício 4 ---")
    
    # Parâmetros ESPECÍFICOS do Exercício 4
    H_k_ex4 = 2000e-6    # 2000 
    beta_ad_ex4 = 0.0    # Sem amortecimento
    p_inlet_dim_ex4 = 0.0 # Pressão nula
    t_max_ex4 = 20.0
    
    # 1. Monta a Rede Hidráulica localmente
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001
    C = AssemblyVectorC(conec, Xno, H_k_ex4, mu) 
    A_net_sparse = sparse.csr_matrix(Assembly(conec, C))
    n_inlet, n_outlet = 0, 5 
    
    # 2. Monta a Membrana localmente
    K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h_ad)
    nm = Nx * Ny
    
    # 3. Acha o 3º Modo
    lambdas, Phi = eigsh(K_mem, k=5, M=M_mem, which='SM')
    idx = np.argsort(lambdas)
    lambdas, Phi = lambdas[idx], Phi[:, idx]
    modo_3 = Phi[:, 2]
    
    fator_Hz = (1.0 / (2.0 * np.pi * R_fisico)) * np.sqrt(sigma / (rho * e_esp))
    freq_analitica_Hz = np.sqrt(lambdas[2]) * fator_Hz
    print(f"-> Frequência Analítica do 3º Modo (Isolada): {freq_analitica_Hz:.2f} Hz")
    
    # 4. Acoplamento (Chama as funções contidas neste mesmo arquivo)
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem, M_mem, A_net_sparse, Nx, Ny, Lx_ad, Ly_ad, R_ad, h_ad, 
        dt, beta_ad_ex4, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet
    )
    
    # 5. Loop no tempo
    p_inlet_adim = p_inlet_dim_ex4 / pref
    w_n = (modo_3 / np.max(np.abs(modo_3))) * 0.1 
    v_n = np.zeros(nm)
    
    tempos = np.arange(0, t_max_ex4 + dt, dt)
    hist_w_centro = []
    no_centro = int((Ny//2) * Nx + (Nx//2))
    
    for t in tempos:
        w_n, v_n, p_n = Resolve_Passo_Tempo(Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, A_net_sparse.shape[0], n_inlet)
        hist_w_centro.append(w_n[no_centro] * (0.01 * R_fisico))
        
    # 6. Processamento dos Picos e Cálculo da Frequência Simulada
    hist_w_array = np.array(hist_w_centro)
    picos, _ = find_peaks(hist_w_array)
    
    if len(picos) > 1:
        tempos_picos = tempos[picos]
        periodos_ad = np.diff(tempos_picos)
        periodo_medio_ad = np.mean(periodos_ad)
        
        freq_simulada_ad = 1.0 / periodo_medio_ad
        freq_simulada_Hz = freq_simulada_ad * fator_Hz
        
        erro_percentual = abs(freq_analitica_Hz - freq_simulada_Hz) / freq_analitica_Hz * 100
        
        print(f"-> Frequência Simulada (Acoplada): {freq_simulada_Hz:.2f} Hz")
        print(f"-> Diferença (Erro devido ao acoplamento): {erro_percentual:.2f}%")
    else:
        print("-> Erro: Não foram encontrados picos suficientes na oscilação.")
        freq_simulada_Hz = 0.0

    # 7. Plotagem dos Resultados do Exercício 4
    plt.figure(figsize=(10, 5))
    plt.plot(tempos, hist_w_centro, 'r-', linewidth=1.5, label='Deslocamento do Nó Central')
    
    if len(picos) > 1:
        plt.plot(tempos[picos], hist_w_array[picos], 'bo', markersize=6, label='Picos Identificados')
        
    plt.title(f"Exercício 4: Decaimento Livre Transiente (3º Modo)\n"
              f"Freq. Analítica: {freq_analitica_Hz:.2f} Hz | Freq. Simulada: {freq_simulada_Hz:.2f} Hz")
    plt.xlabel("Tempo Adimensional")
    plt.ylabel("Deslocamento Central (m)")
    plt.legend(loc="upper right")
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.show()

    

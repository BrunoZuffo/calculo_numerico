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

def calcula_condutancias_quadradas(conec, Xno, H_k, mu):
    conec = np.asarray(conec, dtype=int)
    nc = conec.shape[0]
    C = np.zeros(nc)
    area = H_k ** 2
    D_h = np.sqrt(4.0 * area / np.pi)
    kappa = np.pi * D_h**4 / (128.0 * mu)
    for k in range(nc):
        n1, n2 = conec[k, 0], conec[k, 1]
        x1, y1 = Xno[n1]
        x2, y2 = Xno[n2]
        L_k = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        C[k] = kappa / L_k
    return C

def roda_exercicio_3():
    print("--- Iniciando Exercício 3: Relaxamento Transiente (Corte de Pressão) ---")

    # Parâmetros Fixos do Enunciado
    Nx, Ny = 51, 51
    dt = 0.025
    H_k = 1000e-6          # Canal de 1000 um
    p_inlet_dim = 10000.0  # Pressão inicial de 10 kPa
    
    t_carga = 12.0         # Fase 1: Tempo inflando a membrana
    t_descarga = 18.0      # Fase 2: Tempo após fechar a pressão (relaxamento)
    t_total = t_carga + t_descarga
    
    # Propriedades Físicas
    mu = 5e-4
    mm_to_m = 0.001
    R_fisico = 0.25e-2
    e_esp = 0.1e-3
    sigma = 200.0
    rho = 900.0
    
    # Parâmetros Adimensionais de Referência
    Lx_ad, Ly_ad, R_ad = 2.0, 2.0, 1.0
    sigma_ad, rho_ad, e_ad, beta_ad = 1.0, 1.0, 1.0, 0.1
    h_ad = Lx_ad / (Nx - 1)
    nm = Nx * Ny
    no_centro = (Ny // 2) * Nx + (Nx // 2)

    # Construindo a infraestrutura física (Rede e Membrana)
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * mm_to_m
    n_inlet, n_outlet = 0, 5
    
    C = calcula_condutancias_quadradas(conec, Xno, H_k, mu)
    A_net_sparse = sparse.csr_matrix(Assembly(conec, C))
    np_net = A_net_sparse.shape[0]

    K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h_ad)

    # Montando o acoplamento global
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem, M_mem, A_net_sparse, Nx, Ny, Lx_ad, Ly_ad, R_ad, h_ad, 
        dt, beta_ad, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet
    )

    # Definição de Escalas Físicas para conversão
    wref = 0.01 * R_fisico
    qref = (R_fisico ** 2) * vref

    # Alocação de variáveis para a simulação
    w_n = np.zeros(nm)
    v_n = np.zeros(nm)
    tempos = np.arange(0.0, t_total + 0.5 * dt, dt)
    U_out = U_sparse[n_outlet, :].toarray().ravel()

    # Listas para guardar dados dos gráficos
    hist_t = []
    hist_p_out = []
    hist_q_out = []
    hist_w_centro = []
    hist_volume = []
    hist_pot = []

    # LOOP TEMPORAL
    for t in tempos:
        # AQUI ESTÁ A LÓGICA DO EXERCÍCIO 3:
        # Se o tempo passou de 12, a pressão cai para ZERO instantaneamente.
        if t <= t_carga:
            p_inlet_atual_dim = p_inlet_dim
        else:
            p_inlet_atual_dim = 0.0
            
        p_inlet_adim = p_inlet_atual_dim / pref

        # Resolve o sistema linear acoplado para o instante atual
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, np_net, n_inlet
        )

        # Conversão dos resultados adimensionais para unidades físicas reais
        p_outlet_fisico = p_n[n_outlet] * pref
        q_out_fisico = (h_ad**2 * np.sum(U_out * v_n)) * qref
        w_centro_fisico = w_n[no_centro] * wref
        
        # Cálculo do volume geométrico do reservatório abaulado (m³)
        area_elemento_fisico = (h_ad * R_fisico)**2
        vol_reservatorio = np.sum(w_n * wref) * area_elemento_fisico
        
        # Potência Hidráulica fornecida na entrada (W)
        potencia = p_inlet_atual_dim * q_out_fisico

        # Armazenando nas listas com multiplicadores para melhor escala visual nos eixos
        hist_t.append(t)
        hist_p_out.append(p_outlet_fisico)
        hist_q_out.append(q_out_fisico * 1e9)       # mm³/s
        hist_w_centro.append(w_centro_fisico * 1e6)  # µm
        hist_volume.append(vol_reservatorio * 1e9)   # mm³
        hist_pot.append(potencia * 1e6)              # µW

    # GERANDO OS GRÁFICOS EXIGIDOS
    print("Simulação concluída. Renderizando gráficos...")
    fig, axs = plt.subplots(3, 2, figsize=(14, 10))
    
    # Título principal detalhado conforme a referência
    fig.suptitle(r'Exercício 3 - Rede de canais com $H_k = 1000\ \mu m$ e $p_{inlet} = 10000\ Pa$', fontsize=15)

    # Gráfico 1: Pressão de Saída
    axs[0, 0].plot(hist_t, hist_p_out, 'r-', lw=2)
    axs[0, 0].axvline(x=t_carga, color='k', linestyle='--', alpha=0.5)
    axs[0, 0].set_title(r'Pressão no nó de saída ($p_{outlet}$)')
    axs[0, 0].set_ylabel('Pressão [Pa]')
    axs[0, 0].grid(True, alpha=0.3)

    # Gráfico 2: Vazão de Saída
    axs[0, 1].plot(hist_t, hist_q_out, 'b-', lw=2)
    axs[0, 1].axvline(x=t_carga, color='k', linestyle='--', alpha=0.5)
    axs[0, 1].set_title(r'Vazão no nó de saída ($q_{outlet}$)')
    axs[0, 1].set_ylabel('Vazão [mm³/s]') # Mantido em mm³/s para melhor leitura visual
    axs[0, 1].grid(True, alpha=0.3)

    # Gráfico 3: Deflexão Central
    axs[1, 0].plot(hist_t, hist_w_centro, 'g-', lw=2)
    axs[1, 0].axvline(x=t_carga, color='k', linestyle='--', alpha=0.5)
    axs[1, 0].set_title(r'Deflexão no nó central ($w_0$)')
    axs[1, 0].set_ylabel(r'Deflexão [$\mu$m]')
    axs[1, 0].grid(True, alpha=0.3)

    # Gráfico 4: Volume do Reservatório
    axs[1, 1].plot(hist_t, hist_volume, 'm-', lw=2)
    axs[1, 1].axvline(x=t_carga, color='k', linestyle='--', alpha=0.5)
    axs[1, 1].set_title(r'Volume do reservatório ($V$)')
    axs[1, 1].set_ylabel('Volume [mm³]')
    axs[1, 1].grid(True, alpha=0.3)

    # Gráfico 5: Potência Consumida
    axs[2, 0].plot(hist_t, hist_pot, 'c-', lw=2)
    axs[2, 0].axvline(x=t_carga, color='k', linestyle='--', alpha=0.5)
    axs[2, 0].set_title(r'Potência consumida ($P$)')
    axs[2, 0].set_ylabel(r'Potência [$\mu$W]')
    axs[2, 0].set_xlabel('Tempo adimensional')
    axs[2, 0].grid(True, alpha=0.3)

    # Oculta o sexto quadrante sobressalente
    axs[2, 1].axis('off')

    # Aplicando o eixo X padronizado para os gráficos visíveis
    axs[0, 0].set_xlabel('Tempo adimensional')
    axs[0, 1].set_xlabel('Tempo adimensional')
    axs[1, 0].set_xlabel('Tempo adimensional')
    axs[1, 1].set_xlabel('Tempo adimensional')

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    
    # Salva uma cópia em alta definição para o seu relatório técnico
    plt.savefig('graficos_exercicio3_final.png', dpi=300)
    plt.show()

def roda_exercicio_4(Nx, Ny, h_ad, dt, mu, R_fisico, e_esp, sigma, rho, Lx_ad, Ly_ad, R_ad, sigma_ad, rho_ad, e_ad):
    print("\n--- Iniciando Exercício 4 ---")
    
    H_k_ex4 = 2000e-6    
    beta_ad_ex4 = 0.0    
    p_inlet_dim_ex4 = 0.0 
    t_max_ex4 = 20.0
    
    # 1. Montagem da Rede
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001
    C = AssemblyVectorC(conec, Xno, H_k_ex4, mu) 
    A_net_sparse = sparse.csr_matrix(Assembly(conec, C))
    n_inlet, n_outlet = 0, 5 
    
    # 2. Montagem da Membrana
    K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h_ad)
    nm = Nx * Ny
    
    # 3. Extração do 3º Modo (Garantindo base determinística)
    vetor_inicial = np.ones(nm)
    lambdas, Phi = eigsh(K_mem, k=5, M=M_mem, which='SM', v0=vetor_inicial)

    idx = np.argsort(lambdas)
    lambdas, Phi = lambdas[idx], Phi[:, idx]
    modo_3 = Phi[:, 2]
    
    # Localiza automaticamente o Antinó do modo (ponto de máxima amplitude real)
    no_medicao = np.argmax(np.abs(modo_3))
    
    fator_Hz = (1.0 / (2.0 * np.pi * R_fisico)) * np.sqrt(sigma / (rho * e_esp))
    freq_analitica_Hz = np.sqrt(lambdas[2]) * fator_Hz
    print(f"-> Frequência Analítica do 3º Modo (Isolada): {freq_analitica_Hz:.2f} Hz")
    
    # 4. Acoplamento
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem, M_mem, A_net_sparse, Nx, Ny, Lx_ad, Ly_ad, R_ad, h_ad, 
        dt, beta_ad_ex4, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet
    )
    
    # 5. Evolução Temporal
    p_inlet_adim = p_inlet_dim_ex4 / pref
    w_n = (modo_3 / np.max(np.abs(modo_3))) * 0.1 
    v_n = np.zeros(nm)
    
    tempos = np.arange(0, t_max_ex4 + dt, dt)
    hist_medicao = []
    
    for t in tempos:
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, A_net_sparse.shape[0], n_inlet
        )
        # Extrai deslocamento no antinó
        hist_medicao.append(w_n[no_medicao] * (0.01 * R_fisico))

    # 6. Processamento
    hist_array = np.array(hist_medicao)
    picos, _ = find_peaks(hist_array)
    
    if len(picos) > 1:
        tempos_picos = tempos[picos]
        periodos_ad = np.diff(tempos_picos)
        periodo_medio_ad = np.mean(periodos_ad)
        
        # Frequência cíclica adimensional
        freq_simulada_ad = 1.0 / periodo_medio_ad
        
        # Tempo de referência (t_ref) da adimensionalização
        t_ref = R_fisico * np.sqrt((rho * e_esp) / sigma)
        
        # Conversão correta (sem a divisão de 2*pi atrelada à grandeza angular)
        freq_simulada_Hz = freq_simulada_ad / t_ref
        
        erro_percentual = abs(freq_analitica_Hz - freq_simulada_Hz) / freq_analitica_Hz * 100
        
        print(f"-> Frequência Simulada (Acoplada): {freq_simulada_Hz:.2f} Hz")
        print(f"-> Diferença (Erro Numérico/Acoplamento): {erro_percentual:.2f}%")
    else:
        print("-> Erro: Não foram encontrados picos suficientes na oscilação.")
        freq_simulada_Hz = 0.0

    # 7. Documentação Gráfica
    plt.figure(figsize=(10, 5))
    plt.plot(tempos, hist_medicao, 'r-', linewidth=1.5, label='Deslocamento no Antinó')
    
    if len(picos) > 1:
        plt.plot(tempos[picos], hist_array[picos], 'bo', markersize=6, label='Picos de Oscilação')
        
    plt.title(f"Exercício 4: Decaimento Livre Transiente (3º Modo)\n"
              f"Freq. Analítica: {freq_analitica_Hz:.2f} Hz | Freq. Simulada: {freq_simulada_Hz:.2f} Hz")
    plt.xlabel("Tempo Adimensional")
    plt.ylabel("Deslocamento (m)")
    plt.legend(loc="upper right")
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.tight_layout()
    plt.show()

def roda_exercicio_5(Nx, Ny, h_ad, dt, mu, R_fisico, e_esp, sigma, rho, Lx_ad, Ly_ad, R_ad, sigma_ad, rho_ad, e_ad):
    print("\n--- Iniciando Exercício 5 ---")
    
    H_k_ex5 = 2000e-6
    beta_ad_ex5 = 0.0
    t_max_ex5 = 30.0  
    
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001
    C = AssemblyVectorC(conec, Xno, H_k_ex5, mu)
    A_net_sparse = sparse.csr_matrix(Assembly(conec, C))
    n_inlet, n_outlet = 0, 5 
    
    K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx, Ny, sigma_ad, rho_ad, e_ad, h_ad)
    nm = Nx * Ny
    
    # Extração determinística da frequência alvo
    vetor_inicial = np.ones(nm)
    lambdas, _ = eigsh(K_mem, k=5, M=M_mem, which='SM', v0=vetor_inicial)
    idx = np.argsort(lambdas)
    lambdas = lambdas[idx]
    omega_3_ad = np.sqrt(lambdas[2])
    
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem, M_mem, A_net_sparse, Nx, Ny, Lx_ad, Ly_ad, R_ad, h_ad, 
        dt, beta_ad_ex5, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet
    )
    
    w_n = np.zeros(nm)
    v_n = np.zeros(nm)
    tempos = np.arange(0, t_max_ex5 + dt, dt)
    
    hist_p_out = []
    hist_w_max = []
    
    for t_ad in tempos:
        # Avaliação no instante t_{n+1} para Euler Implícito
        t_next = t_ad + dt
        p_inlet_dim_atual = 5000.0 * np.cos(omega_3_ad * t_next)
        p_inlet_adim = p_inlet_dim_atual / pref
        
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob, M_mem, w_n, v_n, p_inlet_adim, dt, nm, A_net_sparse.shape[0], n_inlet
        )
        
        hist_p_out.append(p_n[n_outlet] * pref)
        # Monitoramento do pico absoluto em vez de um nó isolado
        hist_w_max.append(np.max(np.abs(w_n)) * (0.01 * R_fisico))
        
    # Análise de resultados: Demonstração da não-ressonância assimétrica
    fig, ax1 = plt.subplots(figsize=(10, 5))
    
    ax1.plot(tempos, np.array(hist_w_max) * 1000, 'b-', linewidth=1.5, label='Deslocamento Máximo (mm)')
    ax1.set_xlabel("Tempo Adimensional")
    ax1.set_ylabel("Deslocamento (mm)", color='b')
    ax1.tick_params(axis='y', labelcolor='b')
    ax1.grid(True, linestyle=':', alpha=0.7)
    
    ax2 = ax1.twinx()
    ax2.plot(tempos, hist_p_out, 'r--', alpha=0.6, label='Pressão Outlet (Pa)')
    ax2.set_ylabel("Pressão (Pa)", color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    
    plt.title(f"Exercício 5: Análise de Forçamento na Freq. do 3º Modo ($\omega_3={omega_3_ad:.2f}$ rad/s)")
    fig.tight_layout()
    plt.show()

#Exercício 1 =================================================
def Calcula_Matriz_R(A_net, N1, N2, Lx_ad, Ly_ad, R_ad, h_ad, sigma, rho, e_esp, R_fisico, n_inlet, n_outlet):
#calcula matriz de amortecimento hidraulico
    nm = N1 * N2
    np_net = A_net.shape[0]

    # Adimensionaliza e aplica BC (mesmo procedimento de Monta_Matriz_Global_Acoplada)
    w0   = 0.01 * R_fisico
    pref = sigma * w0 / R_fisico**2
    vref = (w0 / R_fisico) * np.sqrt(sigma / (rho * e_esp))

    A_scaled_lil = (A_net.copy() * (pref / (vref * R_fisico**2))).tolil()
    A_scaled_lil[n_inlet, :] = 0.0
    A_scaled_lil[n_inlet, n_inlet] = 1.0
    A_hat = A_scaled_lil.toarray()   # denso para inversão direta

    # Monta vetores
    xc, yc = Lx_ad / 2.0, Ly_ad / 2.0
    x = np.linspace(0, Lx_ad, N1)
    y = np.linspace(0, Ly_ad, N2)

    uns = np.ones(nm)
    for j in range(N2):
        for i in range(N1):
            dist2 = (x[i] - xc)**2 + (y[j] - yc)**2
            if dist2 >= R_ad**2:
                uns[i + j * N1] = 0.0

    U = np.zeros((np_net, nm))
    U[n_outlet, :] = uns

    # formula R = ĥ² · U^T · A_hat^{-1} · U
    h2  = h_ad**2
    A_inv = np.linalg.inv(A_hat)
    R_mat = h2 * U.T @ A_inv @ U

    return R_mat
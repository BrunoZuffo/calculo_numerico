import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy import sparse
import time

# Importações dos módulos auxiliares
from functionsHT import (ObterCondutividadeFaces_ViaNos, CriarSistemaSolidoCondutividadeVariavel,
                         calcular_temperatura_media_arestas, atualiza_condutancias, 
                         plot_arestas_cromaticas_hidraulics, calcular_viscosidade)

# Importações obrigatórias que estavam faltando
from functions import GeraGrafo, Assembly as AssemblyHidraulico

if __name__ == "__main__":
    
    # =========================================================================
    # 1. MOTOR PRINCIPAL (CONDIÇÕES DE CONTORNO, GEOMETRIA E BASE ISOTÉRMICA)
    # =========================================================================
    print("=" * 80)
    print("1. MOTOR PRINCIPAL: INICIALIZAÇÃO DE PARÂMETROS E REDE DE REFERÊNCIA")
    print("=" * 80)
    
    # Parâmetros geométricos da placa
    Lx, Ly = 0.03, 0.015
    k_0 = 0.25
    fonte_calor = 5.0e5
    TL, TR = 10.0, 30.0
    
    # Parâmetros da inclusão circular
    R_incl = 0.0025
    xincl, yincl = 0.02 + R_incl, 0.0075
    TC = 35.0
    
    # Parâmetros geométricos e operacionais dos microcanais
    Area_canal = 2.5e-7  # 500 um x 500 um
    Q_in_0 = 1.0e-7      # 0.1 mL/s (Inlet nó 0)
    Q_in_175 = 1.0e-6    # 1.0 mL/s (Inlet nó 175)
    
    # Construção da topologia da rede hidráulica fractal
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001       # Conversão de mm para metros
    Xno[:, 1] += 0.5 * Ly   # Centralização geométrica no eixo Y da placa
    
    # Nós de controle de contorno hidráulico
    n_inlet_0 = 0
    n_inlet_175 = 175
    n_outlet = 5 
    
    # Resolução da hidráulica nominal invariante (Caso Isotérmico de Referência a T = 20°C)
    mu_iso = calcular_viscosidade(20.0)
    D_h = np.sqrt(4 * Area_canal / np.pi)
    
    C_iso = np.zeros(len(conec))
    for k in range(len(conec)):
        Lk = np.linalg.norm(Xno[int(conec[k,1])] - Xno[int(conec[k,0])])
        C_iso[k] = (np.pi * D_h**4) / (128 * mu_iso * Lk)
        
    A_h = AssemblyHidraulico(conec, C_iso)
    A_h_tilde = A_h.copy()
    A_h_tilde[n_outlet, :] = 0.0
    A_h_tilde[n_outlet, n_outlet] = 1.0
    
    b_h = np.zeros(len(Xno))
    b_h[n_inlet_0] = Q_in_0
    b_h[n_inlet_175] = Q_in_175
    
    pressures_iso = np.linalg.solve(A_h_tilde, b_h)
    P_max_iso = np.max(pressures_iso)
    print(f"Pressão Máxima Isotérmica Base calculada: {P_max_iso:.2f} Pa\n")


    # =========================================================================
    # 2. EXERCÍCIOS DA PARTE 1 (SEÇÃO 4.2.1 - ACOPLAMENTO HIDRO-TÉRMICO)
    # =========================================================================
    print("=" * 80)
    print("2. EXERCÍCIOS DA PARTE 1 (SEÇÃO 4.2.1): COUPLING E VISCOSIDADE VARIÁVEL")
    print("=" * 80)
    
    # 2.1 Resolução do campo de temperaturas nominal da placa para amostragem (Malha Refinada)
    Nx_ref, Ny_ref = 241, 121
    d_max_ref = 0.0005
    
    x_ref_coords = np.linspace(0, Lx, Nx_ref)
    TB_ref = 10.0 + 20.0 * (x_ref_coords / Lx)
    TT_ref = 10.0 + 20.0 * (x_ref_coords / Lx)
    
    k_e_ref, k_n_ref = ObterCondutividadeFaces_ViaNos(Nx_ref, Ny_ref, Lx, Ly, Xno, conec, d_max_ref, k_0)
    A_s_ref, b_s_ref = CriarSistemaSolidoCondutividadeVariavel(
        Nx_ref, Ny_ref, Lx, Ly, k_e_ref, k_n_ref, TL, TR, TB_ref, TT_ref, fonte_calor, R_incl, xincl, yincl, TC
    )
    T_solid_ref = spsolve(sparse.csr_matrix(A_s_ref), b_s_ref)
    
    # Tarefa 4.2.1.2 & 4.2.1.3: Teste de convergência da quadratura (N = 10, 100, 1000)
    print("Executando amostragem ao longo das arestas via RegularGridInterpolator:")
    subintervalos = [10, 100, 1000]
    T_arestas_final = None
    
    for N in subintervalos:
        T_arestas, t_exec = calcular_temperatura_media_arestas(
            conec, Xno, T_solid_ref, Nx_ref, Ny_ref, Lx, Ly, num_subintervalos=N
        )
        if N == 1000:
            T_arestas_final = T_arestas
        print(f" -> N = {N:<4} | Tempo de execução: {t_exec:.4f} s | Temperatura Média Global: {np.mean(T_arestas):.4f} °C")
        
    # Tarefa 4.2.1.4: Atualização das propriedades moleculares e recálculo das pressões
    C_acoplado = atualiza_condutancias(conec, Xno, T_arestas_final, Area_canal)
    
    A_h_mod = AssemblyHidraulico(conec, C_acoplado)
    A_h_mod_tilde = A_h_mod.copy()
    A_h_mod_tilde[n_outlet, :] = 0.0
    A_h_mod_tilde[n_outlet, n_outlet] = 1.0
    
    pressures_mod = np.linalg.solve(A_h_mod_tilde, b_h)
    P_max_mod = np.max(pressures_mod)
    
    print("-" * 80)
    print(f"Pressão Máxima (Isotérmica):          {P_max_iso:.2f} Pa")
    print(f"Pressão Máxima (Acoplada Corrigida):  {P_max_mod:.2f} Pa")
    print(f"Variação da Perda de Carga Realizada: {((P_max_mod - P_max_iso)/P_max_iso)*100:.2f}%")
    print("-" * 80)
    
    # Renderização dos perfis de controle térmico do fluido
    plot_arestas_cromaticas_hidraulics(conec, Xno, T_arestas_final, titulo="Distribuição Térmica Integrada nas Arestas da Rede (N=1000)")
    print("\n")


    # =========================================================================
    # 3. EXERCÍCIOS DA PARTE 2 (SEÇÃO 4.3.3 - ANÁLISE PARAMÉTRICA DO SÓLIDO)
    # =========================================================================
    print("=" * 80)
    print("3. EXERCÍCIOS DA PARTE 2 (SEÇÃO 4.3.3): ANÁLISE PARAMÉTRICA")
    print("=" * 80)
    
    # --- TAREFA 4.3.3.1: Estudo de Convergência de Malha e d_max ---
    print("Iniciando varredura paramétrica de malhas e distâncias de corte:")
    malhas = [(61, 31), (121, 61), (241, 121)]
    d_max_valores = [0.00025, 0.0005, 0.0010]
    
    print(f"{'Malha (Nx x Ny)':<18} | {'d_max (mm)':<12} | {'T_max (°C)':<12} | {'Tempo K (s)':<12} | {'Tempo Solv (s)':<12}")
    print("-" * 80)
    
    for (Nx, Ny) in malhas:
        x_coords = np.linspace(0, Lx, Nx)
        TB = 10.0 + 20.0 * (x_coords / Lx)
        TT = 10.0 + 20.0 * (x_coords / Lx)
        
        for d_max_val in d_max_valores:
            t_k_ini = time.perf_counter()
            
            # CORREÇÃO: Utilizar Nx, Ny e d_max_val dinâmicos do laço, e não os de referência (_ref)
            k_e, k_n = ObterCondutividadeFaces_ViaNos(Nx, Ny, Lx, Ly, Xno, conec, d_max_val, k_0)
            
            t_k_fim = time.perf_counter()
            
            t_s_ini = time.perf_counter()
            
            # CORREÇÃO: Passar k_e e k_n diretamente
            A_s, b_s = CriarSistemaSolidoCondutividadeVariavel(
                Nx, Ny, Lx, Ly, k_e, k_n, TL, TR, TB, TT, fonte_calor, R_incl, xincl, yincl, TC
            )
            T_solid_flat = spsolve(sparse.csr_matrix(A_s), b_s)
            t_s_fim = time.perf_counter()
            
            print(f"{f'{Nx} x {Ny}':<18} | {d_max_val*1000:<12.3f} | {np.max(T_solid_flat):<12.2f} | {t_k_fim - t_k_ini:<12.4f} | {t_s_fim - t_s_ini:<12.4f}")
            
    print("-" * 80)
    
    # --- ESPAÇO PARA A TAREFA 4.3.3.2: EFEITO DO TERMO FONTE PERTURBADO (S_p) ---
    # TODO: Implementar a modificação do vetor de forças 'b_s' adicionando a contribuição
    # das intensidades I_j das correntes associadas às ramificações (homogênea vs heterogênea).
    print("\n[INFO] Espaço reservado para a implementação da Tarefa 4.3.3.2 (Termo Fonte Epigenético S_p).")
    # =========================================================================
    # IMPLEMENTAR AQUI O LAÇO PARA S_0 = [+/-1e5, +/-5e5, +/-1e6] E VARIAÇÃO DE I_j
    # =========================================================================
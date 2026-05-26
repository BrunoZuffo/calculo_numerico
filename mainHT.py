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
    print("\n" + "="*80)
    print("3. EXERCÍCIOS DA PARTE 2 (SEÇÃO 4.3.3): ANÁLISE PARAMÉTRICA")
    print("="*80)

    # --- TAREFA 4.3.3.1: Estudo de Convergência de Malha e d_max ---
    print("\n" + "="*50)
    print("INICIANDO EXERCÍCIO 1: ANÁLISE PARAMÉTRICA DE D_MAX E MALHAS")
    print("="*50)
    
    lista_malhas = [(61, 31), (121, 61), (241, 121)]
    lista_dmax = [0.00025, 0.0005, 0.001]

    for (Nx, Ny) in lista_malhas:
        for d_max_val in lista_dmax:
            print(f"\n-> Analisando Malha ({Nx}x{Ny}) com d_max = {d_max_val}")
            t_inicio = time.time()
            
            x_coords = np.linspace(0, Lx, Nx)
            y_coords = np.linspace(0, Ly, Ny)
            TB = 10.0 + 20.0 * (x_coords / Lx)
            TT = 10.0 + 20.0 * (x_coords / Lx)
            
            # Etapa 1: Calcula o campo de condutividade nas faces (FUNÇÃO CORRIGIDA)
            t0 = time.time()
            k_e, k_n = ObterCondutividadeFaces_ViaNos(Nx, Ny, Lx, Ly, Xno, conec, d_max_val, k_0)
            t_k = time.time() - t0
            
            # Etapa 2: Montagem do sistema
            t0 = time.time()
            A_s, b_s = CriarSistemaSolidoCondutividadeVariavel(Nx, Ny, Lx, Ly, k_e, k_n, TL, TR, TB, TT, fonte_calor, R_incl, xincl, yincl, TC)
            A_s_sparse = sparse.csr_matrix(A_s)
            t_assembly = time.time() - t0
            
            # Etapa 3: Resolução do sistema
            t0 = time.time()
            T_solid_flat = spsolve(A_s_sparse, b_s)
            t_solve = time.time() - t0
            
            t_total = time.time() - t_inicio
            T_max = np.max(T_solid_flat)
            T_placa_2D = T_solid_flat.reshape((Ny, Nx))
            
            print(f"   [Tempo] Calc k_faces: {t_k:.3f}s | Assembly: {t_assembly:.3f}s | Solve: {t_solve:.3f}s | Total: {t_total:.3f}s")
            print(f"   [Resultado] Temperatura Máxima na Placa: {T_max:.2f} °C")
            
            # --- GERAÇÃO DOS GRÁFICOS ---
            fig = plt.figure(figsize=(15, 4))
            
            # 1. Mapa de contornos 2D
            ax1 = plt.subplot(1, 3, 1)
            cf = ax1.contourf(x_coords * 1000, y_coords * 1000, T_placa_2D, levels=50, cmap='jet')
            plt.colorbar(cf, ax=ax1, label='T (°C)')
            ax1.set_title(f'Mapa 2D ($d_{{max}}$={d_max_val})')
            ax1.set_xlabel('X (mm)')
            ax1.set_ylabel('Y (mm)')
            
            # 2. Perfil Horizontal (Eixo Central em y = Ly/2)
            ax2 = plt.subplot(1, 3, 2)
            idx_y_mid = Ny // 2
            ax2.plot(x_coords * 1000, T_placa_2D[idx_y_mid, :], 'r-', label=f'Y = {y_coords[idx_y_mid]*1000:.1f} mm')
            ax2.set_title('Perfil Horizontal')
            ax2.set_xlabel('X (mm)')
            ax2.set_ylabel('T (°C)')
            ax2.grid(True)
            ax2.legend()
            
            # 3. Perfil Vertical (Eixo Central em x = Lx/2)
            ax3 = plt.subplot(1, 3, 3)
            idx_x_mid = Nx // 2
            ax3.plot(y_coords * 1000, T_placa_2D[:, idx_x_mid], 'b-', label=f'X = {x_coords[idx_x_mid]*1000:.1f} mm')
            ax3.set_title('Perfil Vertical')
            ax3.set_xlabel('Y (mm)')
            ax3.set_ylabel('T (°C)')
            ax3.grid(True)
            ax3.legend()
            
            plt.tight_layout()
            plt.show()
    
    # --- ESPAÇO PARA A TAREFA 4.3.3.2: EFEITO DO TERMO FONTE PERTURBADO (S_p) ---
    # TODO: Implementar a modificação do vetor de forças 'b_s' adicionando a contribuição
    # das intensidades I_j das correntes associadas às ramificações (homogênea vs heterogênea).
    print("\n[INFO] Espaço reservado para a implementação da Tarefa 4.3.3.2 (Termo Fonte Epigenético S_p).")
    # =========================================================================
    # IMPLEMENTAR AQUI O LAÇO PARA S_0 = [+/-1e5, +/-5e5, +/-1e6] E VARIAÇÃO DE I_j
    # =========================================================================
import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse.linalg import spsolve
from scipy import sparse
from scipy.spatial import cKDTree
import time
from functionsHT import calcular_temperatura_media_arestas, plot_arestas_cromaticas_hidraulics, CalcularK_Face, ObterCondutividadeFaces, CriarSistemaSolidoCondutividadeVariavel

# Importações dos teus módulos existentes
from functions import GeraGrafo, AssemblyVectorC, Assembly as AssemblyHidraulico, calc_vazao

# =====================================================================
# 1. FUNÇÕES DO PROFESSOR (mapea_pontos_prox.py integradas)
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
# 2. FUNÇÕES DE ACOPLAMENTO TÉRMICO E ESTRUTURAÇÃO DA PLACA
# =====================================================================
def CriarSistemaSolidoBase(Nx, Ny, Lx, Ly, k_solid, TL, TR, TB, TT, fonte_calor):
    """Monta a matriz de condução pura do sólido (Diferenças Finitas)"""
    nunk = Nx * Ny
    A = np.zeros((nunk, nunk))
    b = np.zeros(nunk)
    hx = Lx / (Nx - 1)
    hy = Ly / (Ny - 1)
    
    for i in range(Nx):
        for j in range(Ny):
            Ic = i + j * Nx
            
            # Condições de contorno de Dirichlet nas bordas da placa
            if i == 0:
                A[Ic, Ic] = 1.0
                b[Ic] = TL
            elif i == Nx - 1:
                A[Ic, Ic] = 1.0
                b[Ic] = TR
            elif j == 0:
                A[Ic, Ic] = 1.0
                b[Ic] = TB[i]
            elif j == Ny - 1:
                A[Ic, Ic] = 1.0
                b[Ic] = TT[i]
            else:
                # Nós internos: -k * Laplacian(T) = fonte_calor
                Ie, Iw = (i + 1) + j * Nx, (i - 1) + j * Nx
                In, Is = i + (j + 1) * Nx, i + (j - 1) * Nx
                
                A[Ic, Ic] = 2 * k_solid / hx**2 + 2 * k_solid / hy**2
                A[Ic, Ie] = -k_solid / hx**2
                A[Ic, Iw] = -k_solid / hx**2
                A[Ic, In] = -k_solid / hy**2
                A[Ic, Is] = -k_solid / hy**2
                b[Ic] = fonte_calor
                
    return A, b

def CalcularComprimentos(Xno, conec):
    Nc = len(conec)
    L_canos = np.zeros(Nc)
    for k, (n1, n2) in enumerate(conec):
        dx = Xno[n2, 0] - Xno[n1, 0]
        dy = Xno[n2, 1] - Xno[n1, 1]
        L_canos[k] = np.sqrt(dx**2 + dy**2)
    return L_canos

def ObterTemperaturaSolidoNosCanais(conec, T_solid_flat, mapa_proximidade):
    """Calcula a temperatura média do sólido que envolve cada canal"""
    Nc = len(conec)
    T_s_canais = np.zeros(Nc)
    contagem_nos_por_canal = np.zeros(Nc)
    
    for id_no, arestas_proximas in enumerate(mapa_proximidade):
        for id_aresta, dist in arestas_proximas:
            T_s_canais[id_aresta] += T_solid_flat[id_no]
            contagem_nos_por_canal[id_aresta] += 1
            
    for k in range(Nc):
        if contagem_nos_por_canal[k] > 0:
            T_s_canais[k] /= contagem_nos_por_canal[k]
        else:
            T_s_canais[k] = np.mean(T_solid_flat)
            
    return T_s_canais, contagem_nos_por_canal

def AssemblyTermicoFluido(Xno, conec, Q, L_canos, P_canos, T_s_canais, rho, cp, hc, T_inlet, n_inlet):
    """Monta o sistema térmico para a temperatura do fluido nos canais (Upwind)"""
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
            Af[n_out, n_out] += rho * cp * abs_qk + hc * Pk * Lk
            Af[n_out, n_in]  += -rho * cp * abs_qk
            bf[n_out]        += hc * Pk * Lk * T_s_canais[k]
            
    Af[n_inlet, :] = 0.0
    Af[n_inlet, n_inlet] = 1.0
    bf[n_inlet] = T_inlet
    return Af, bf

def ModificarSistemaSolido(A_s_base, b_s_base, conec, Q, Tf, L_canos, P_canos, hc, hx, hy, Nx, Ny, mapa_proximidade, contagem_nos_por_canal):
    """Adiciona o termo de troca convectiva volumétrica nos nós internos do sólido"""
    A_s = A_s_base.copy().tolil() if sparse.issparse(A_s_base) else A_s_base.copy()
    b_s = b_s_base.copy()
    
    for id_no, arestas_proximas in enumerate(mapa_proximidade):
        i = id_no % Nx
        j = id_no // Nx
        
        # Mantém intactas as condições de contorno originais da placa
        if i == 0 or i == Nx-1 or j == 0 or j == Ny-1:
            continue
            
        for id_aresta, dist in arestas_proximas:
            qk = Q[id_aresta]
            n1, n2 = conec[id_aresta]
            n_in = n1 if qk >= 0 else n2
            
            Tf_canal = Tf[n_in]
            Lk = L_canos[id_aresta]
            Pk = P_canos[id_aresta]
            num_nos_mapeados = contagem_nos_por_canal[id_aresta]
            
            if num_nos_mapeados > 0:
                # O calor trocado é distribuído proporcionalmente ao volume da célula (hx * hy)
                fator_volumetrico = (hc * Pk * Lk) / (num_nos_mapeados * hx * hy)
                A_s[id_no, id_no] += fator_volumetrico
                b_s[id_no]        += fator_volumetrico * Tf_canal
                
    return A_s.tocsc() if sparse.issparse(A_s) else A_s

# =====================================================================
# 3. CONFIGURAÇÃO E EXECUÇÃO PRINCIPAL
# =====================================================================
if __name__ == "__main__":
    # Parâmetros Geométricos e Térmicos da Placa
    Nx, Ny = 51, 26
    Lx, Ly = 0.02, 0.01  # metros (2cm x 1cm)
    hx, hy = Lx / (Nx - 1), Ly / (Ny - 1)
    
    k_solid = 0.25       
    rho, cp, hc = 1000.0, 4184.0, 3000.0
    T_inlet = 20.0
    TL, TR = 10.0, 30.0
    fonte_calor = 5.0e5

    # 3.1 Resolução Hidráulica da Rede
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001  # Conversão da escala da malha para metros
    n_inlet, n_outlet = 0, len(Xno) - 1

    # Correção do vetor C para garantir que seja unidimensional puro
    C = AssemblyVectorC(conec, Xno)
    C = np.array(C).flatten()  

    A_h = AssemblyHidraulico(conec, C)
    A_h_tilde = A_h.copy()
    A_h_tilde[n_outlet, :] = 0.0
    A_h_tilde[n_outlet, n_outlet] = 1.0

    b_h = np.zeros(len(Xno))
    b_h[n_inlet] = 1.0e-7  # Vazão imposta na entrada

    pressures = np.linalg.solve(A_h_tilde, b_h)
    
    # CORREÇÃO DA ORDEM: conec, C, pressures conforme a assinatura do teu arquivo
    Q = calc_vazao(conec, C, pressures)  

    L_canos = CalcularComprimentos(Xno, conec)
    D_canos = np.ones(len(conec)) * 100e-6  # canais com diâmetro de 100 µm
    P_canos = np.pi * D_canos

    # 3.2 Geração do Mapa de Proximidade (Algoritmo do Professor)
    d_max = 1.5 * np.sqrt(hx**2 + hy**2) 
    print(f"Gerando mapa de proximidade via cKDTree (d_max = {d_max:.4e} m)...")
    mapa_proximidade = CreateMapDistance(Lx, Ly, Nx, Ny, Xno, conec, d_max)

    # 3.3 Construção da Matriz de Condução Base da Placa
    x_coords = np.linspace(0, Lx, Nx)
    TB = 10 + 20 * (x_coords / Lx)
    TT = 10 + 20 * (x_coords / Lx)
    A_s_base, b_s_base = CriarSistemaSolidoBase(Nx, Ny, Lx, Ly, k_solid, TL, TR, TB, TT, fonte_calor)
    A_s_base_sparse = sparse.csc_matrix(A_s_base)

    # 3.4 Loop Iterativo de Acoplamento Térmico (Gauss-Seidel entre Sistemas)
    nunk = Nx * Ny
    T_solid_flat = np.ones(nunk) * T_inlet 
    T_fluid = np.ones(len(Xno)) * T_inlet
    TOL, MAX_ITER = 1e-5, 200
    erro, iteracao = 1.0, 0

    print("\nIniciando acoplamento térmico...")
    while erro > TOL and iteracao < MAX_ITER:
        T_solid_old = T_solid_flat.copy()
        
        # Passo A: Fluido vê o sólido através do mapa e calcula suas novas temperaturas
        T_s_canais, contagem_nos = ObterTemperaturaSolidoNosCanais(conec, T_solid_flat, mapa_proximidade)
        Af, bf = AssemblyTermicoFluido(Xno, conec, Q, L_canos, P_canos, T_s_canais, rho, cp, hc, T_inlet, n_inlet)
        T_fluid = np.linalg.solve(Af, bf)
        
        # Passo B: Sólido recebe a energia dos canais e calcula sua nova distribuição térmica
        A_s_mod = ModificarSistemaSolido(A_s_base_sparse, b_s_base, conec, Q, T_fluid, L_canos, P_canos, hc, hx, hy, Nx, Ny, mapa_proximidade, contagem_nos)
        T_solid_flat = spsolve(A_s_mod, b_s_base)
        
        erro = np.max(np.abs(T_solid_flat - T_solid_old))
        iteracao += 1
        if iteracao % 5 == 0 or erro <= TOL:
            print(f"Iteração: {iteracao:<4} | Erro Máximo: {erro:.4e}")

    # 3.5 Visualização do Resultado Final
    T_placa_2D = T_solid_flat.reshape((Ny, Nx))
    plt.figure(figsize=(10, 5))
    plt.contourf(x_coords * 1000, np.linspace(0, Ly, Ny) * 1000, T_placa_2D, levels=50, cmap='jet')
    plt.colorbar(label='Temperatura (°C)')
    plt.title('Distribuição Térmica Final na Placa (Com Mapeamento por KD-Tree)')
    plt.xlabel('X (mm)')
    plt.ylabel('Y (mm)')
    plt.tight_layout()
    plt.show()

# Supondo que 'conec', 'Xno' e 'T_fluid' (temperatura nos nós do fluido) já existam no seu main...
cenarios = [1, 10, 100, 1000]
resultados_T = {}

print("\n--- ANÁLISE DE CONVERGÊNCIA DE QUADRATURA (CAP 4, EX 3) ---")
for N in cenarios:
    T_medias, tempo = calcular_temperatura_media_arestas(conec, Xno, T_fluid, num_subintervalos=N)
    resultados_T[N] = T_medias
    print(f"Subintervalos: {N:<4} | Tempo: {tempo:.6f} s | Média Global: {np.mean(T_medias):.4f} °C")

# Estratégia matemática para estimativa do erro (Critério de Cauchy contra N=1000)
print("\n--- ESTIMATIVA DE ERRO ABSOLUTO MÉDIO ---")
for N in [1, 10, 100]:
    erro_estimado = np.mean(np.abs(resultados_T[N] - resultados_T[1000]))
    print(f"Configuração N={N:<3} -> Erro Médio Estimado: {erro_estimado:.6e} °C")

# Plota o resultado cromático final para a melhor configuração (N=1000)
plot_arestas_cromaticas_hidraulics(
    conec, Xno, 
    T_nodes=T_fluid, 
    T_arestas=resultados_T[1000], 
    titulo="Temperatura Média nas Arestas (Trapézio Composto N=1000)"
)

# EXERCICIO 1 PARTE 2

print("\n" + "="*50)
print("INICIANDO EXERCÍCIO 1: ANÁLISE PARAMÉTRICA DE D_MAX E MALHAS")
print("="*50)
# Parâmetros atualizados da Seção 4.1
Lx, Ly = 0.03, 0.015         # 3cm x 1.5cm
k_0 = 0.25                   
fonte_calor = 5.0e5          
TL, TR = 10.0, 30.0
R_incl = 0.0025              # Raio de 0.25 cm
xincl, yincl = 0.02 + R_incl, 0.0075  # Centro em (2+R, 0.75) cm
TC = 35.0                    # Temperatura do círculo
# Cenários exigidos
lista_malhas = [(61, 31), (121, 61), (241, 121)]
lista_dmax = [0.00025, 0.0005, 0.001]
# Obtenção da Rede (assumimos que foi gerada na Seção 3.1 com levels=3)
# Reescale a rede se necessário para o novo Lx, Ly da placa
Xno, conec = GeraGrafo(levels=3)
Xno = Xno * 0.001
Xno[:,1] += 0.5 * Ly  # Centralizar em Y se necessário (depende da sua rede)
for (Nx, Ny) in lista_malhas:
    for d_max_val in lista_dmax:
        print(f"\n-> Analisando Malha ({Nx}x{Ny}) com d_max = {d_max_val}")
        t_inicio = time.time()
        
        # Borda Superior e Inferior dinâmicas: T = 10 + 20*(x/Lx)
        x_coords = np.linspace(0, Lx, Nx)
        y_coords = np.linspace(0, Ly, Ny)
        TB = 10.0 + 20.0 * (x_coords / Lx)
        TT = 10.0 + 20.0 * (x_coords / Lx)
        # Etapa 1: Calcula o campo de condutividade nas faces
        t0 = time.time()
        k_e, k_n = ObterCondutividadeFaces(Nx, Ny, Lx, Ly, Xno, conec, d_max_val, k_0)
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





# =====================================================================
    # EXERCÍCIO 1 - SEÇÃO 4.3.3: INFLUÊNCIA DA CONDUTIVIDADE VARIÁVEL
    # =====================================================================
    import time
    from functionsHT import ObterCondutividadeFaces, CriarSistemaSolidoCondutividadeVariavel

    print("\n" + "="*50)
    print("INICIANDO EXERCÍCIO 1: ANÁLISE PARAMÉTRICA DE D_MAX E MALHAS")
    print("="*50)

    # Parâmetros atualizados da Seção 4.1
    Lx, Ly = 0.03, 0.015         # 3cm x 1.5cm
    k_0 = 0.25                   
    fonte_calor = 5.0e5          
    TL, TR = 10.0, 30.0
    R_incl = 0.0025              # Raio de 0.25 cm
    xincl, yincl = 0.02 + R_incl, 0.0075  # Centro em (2+R, 0.75) cm
    TC = 35.0                    # Temperatura do círculo

    # Cenários exigidos
    lista_malhas = [(61, 31), (121, 61), (241, 121)]
    lista_dmax = [0.00025, 0.0005, 0.001]

    # Obtenção da Rede (assumimos que foi gerada na Seção 3.1 com levels=3)
    # Reescale a rede se necessário para o novo Lx, Ly da placa
    Xno, conec = GeraGrafo(levels=3)
    Xno = Xno * 0.001
    Xno[:,1] += 0.5 * Ly  # Centralizar em Y se necessário (depende da sua rede)

    for (Nx, Ny) in lista_malhas:
        for d_max_val in lista_dmax:
            print(f"\n-> Analisando Malha ({Nx}x{Ny}) com d_max = {d_max_val}")
            t_inicio = time.time()
            
            # Borda Superior e Inferior dinâmicas: T = 10 + 20*(x/Lx)
            x_coords = np.linspace(0, Lx, Nx)
            y_coords = np.linspace(0, Ly, Ny)
            TB = 10.0 + 20.0 * (x_coords / Lx)
            TT = 10.0 + 20.0 * (x_coords / Lx)

            # Etapa 1: Calcula o campo de condutividade nas faces
            t0 = time.time()
            k_e, k_n = ObterCondutividadeFaces(Nx, Ny, Lx, Ly, Xno, conec, d_max_val, k_0)
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

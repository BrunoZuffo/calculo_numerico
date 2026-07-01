import time
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from functions import GeraGrafo, Assembly, calc_vazao, calc_potencia
from functionsM import BuildMatrizes_Eigen_Circular
from functionsHT import (
    ObterCondutividadeFaces_ViaNos,
    CriarSistemaSolidoCondutividadeVariavel,
    calcular_temperatura_media_arestas,
    calcular_viscosidade,
    atualiza_condutancias,
)
from functionsHM import (
    Monta_Matriz_Global_Acoplada,
    Resolve_Passo_Tempo,
)

_integral_trapezoidal = getattr(np, "trapezoid", None) or np.trapz


# =======================================================================================
# CÓDIGO BASE - GÊMEO DIGITAL
# =======================================================================================

# ========================================================================
# 1. SUBSISTEMA 1 - PLACA TERMICA
#    (condutividade efetiva modificada pela proximidade dos microcanais)
# ========================================================================
def ConstroiPlacaTermica(Xno, conec, Lx, Ly, Nx, Ny, k0, fonte_calor,
                          TL, TR, R_incl, xincl, yincl, TC, d_max):
    """
    Monta e resolve o campo de temperatura na placa solida, considerando
    a condutividade termica efetiva modificada pela proximidade dos
    microcanais da rede hidraulica (Equacao 4.3 do PDF).

    Condicoes de contorno (Secao 6.2):
        TL = 10 C (lado esquerdo), TR = 30 C (lado direito)
        TB = TT = 10 + 20*(x/Lx)   (base e topo)
        Regiao circular de raio R_incl, centrada em (xincl, yincl),
        com temperatura fixada em TC.

    Retorna o campo de temperatura achatado (T_solid_flat, ordenado por
    Ic = i + j*Nx) e as matrizes de condutividade efetiva nas faces
    leste (k_e) e norte (k_n).
    """
    x_coords = np.linspace(0.0, Lx, Nx)

    # Condicao de contorno de temperatura na base e no topo (Secao 6.2)
    TB = 10.0 + 20.0 * (x_coords / Lx)
    TT = 10.0 + 20.0 * (x_coords / Lx)

    # Condutividade efetiva nas faces, modificada pela proximidade dos canais
    k_e, k_n = ObterCondutividadeFaces_ViaNos(Nx, Ny, Lx, Ly, Xno, conec, d_max, k0)

    # Montagem e resolucao do sistema linear de conducao de calor
    A_solid, b_solid = CriarSistemaSolidoCondutividadeVariavel(
        Nx, Ny, Lx, Ly, k_e, k_n, TL, TR, TB, TT, fonte_calor,
        R_incl, xincl, yincl, TC
    )
    T_solid_flat = np.linalg.solve(A_solid, b_solid)

    return T_solid_flat, k_e, k_n


# ========================================================================
# 2. SUBSISTEMA 2 - REDE HIDRAULICA
#    (condutancias atualizadas pela temperatura local de cada canal)
# ========================================================================
def AtualizaRedeViaTemperatura(conec, Xno, T_solid_flat, Nx, Ny, Lx, Ly, H_k,
                                num_subintervalos=100):
    """
    A partir do campo de temperatura da placa, calcula:
        (a) a temperatura media de cada canal, via integral de linha
            do campo interpolado ao longo de cada aresta da rede;
        (b) a viscosidade dinamica local do fluido, segundo a relacao
            empirica mu(T) da Secao 6.2;
        (c) a condutancia hidraulica atualizada de cada canal, para
            secao transversal quadrada de lado H_k (area = H_k**2).
    """
    Area_canal = H_k ** 2

    T_arestas, _tempo_calc = calcular_temperatura_media_arestas(
        conec, Xno, T_solid_flat, Nx, Ny, Lx, Ly,
        num_subintervalos=num_subintervalos
    )

    mu_arestas = calcular_viscosidade(T_arestas)

    C = atualiza_condutancias(conec, Xno, T_arestas, Area_canal)

    return C, T_arestas, mu_arestas


# ========================================================================
# 3. SUBSISTEMA 3 - MEMBRANA ELASTICA
# ========================================================================
def ConstroiMembrana(Nx_m, Ny_m, Lx_ad, sigma_ad, rho_ad, e_ad):
    """
    Monta as matrizes de rigidez (K) e massa (M) da membrana elastica
    circular adimensionalizada, sobre uma malha Nx_m x Ny_m no dominio
    quadrado [0, Lx_ad] x [0, Ly_ad] (Capitulo 3).
    """
    h_ad = Lx_ad / (Nx_m - 1)
    K_mem, M_mem = BuildMatrizes_Eigen_Circular(Nx_m, Ny_m, sigma_ad, rho_ad, e_ad, h_ad)
    return K_mem, M_mem, h_ad


# ========================================================================
# 4. MONTAGEM COMPLETA DO GEMEO DIGITAL
# ========================================================================

def MontaGemeoDigital(params):
    """
    Executa o pipeline completo de montagem do Gemeo Digital (GD),
    obedecendo a dependencia sequencial estabelecida na Secao 6.1:

        Placa Termica -> Viscosidade/Condutancia dos Canais ->
        Rede Hidraulica -> Acoplamento Monolitico com a Membrana

    `params` e um dicionario com todos os parametros fisicos e
    numericos nominais (ver mainGD.py, Secao 1).

    Retorna um dicionario `GD` contendo todos os objetos montados,
    pronto para ser reaproveitado nos exercicios de investigacao das
    Secoes 6.3 e 6.4.
    """
    # ---- 1. Topologia da rede hidraulica -------------------------------
    Xno, conec = GeraGrafo(levels=params['complex_level'])
    Xno = Xno * params['mm_to_m']
    # Centraliza a rede no eixo Y da placa termica (mesma convencao do
    # acoplamento Hidraulico-Termico do Capitulo 4).
    Xno[:, 1] += 0.5 * params['Ly']

    n_inlet = params['n_inlet']
    n_outlet = params['n_outlet']

    # ---- 2. Placa termica (condutividade modificada pelos canais) -----
    T_solid_flat, k_e, k_n = ConstroiPlacaTermica(
        Xno, conec, params['Lx'], params['Ly'], params['Nx_t'], params['Ny_t'],
        params['k0'], params['fonte_calor'], params['TL'], params['TR'],
        params['R_incl'], params['xincl'], params['yincl'], params['TC'],
        params['d_max']
    )

    # ---- 3. Atualizacao termica das condutancias da rede ---------------
    C, T_arestas, mu_arestas = AtualizaRedeViaTemperatura(
        conec, Xno, T_solid_flat, params['Nx_t'], params['Ny_t'],
        params['Lx'], params['Ly'], params['H_k']
    )

    A_net = Assembly(conec, C)
    A_net_sparse = sparse.csr_matrix(A_net)
    np_net = A_net_sparse.shape[0]

    # ---- 4. Membrana elastica ------------------------------------------
    K_mem, M_mem, h_ad = ConstroiMembrana(
        params['Nx_m'], params['Ny_m'], params['Lx_ad'],
        params['sigma_ad'], params['rho_ad'], params['e_ad']
    )
    nm = params['Nx_m'] * params['Ny_m']

    # ---- 5. Acoplamento monolitico Hidraulico-Mecanico ------------------
    Aglob, U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        K_mem, M_mem, A_net_sparse,
        params['Nx_m'], params['Ny_m'], params['Lx_ad'], params['Ly_ad'],
        params['R_ad'], h_ad, params['dt'], params['beta_ad'],
        params['sigma'], params['rho'], params['e_esp'], params['R_fisico'],
        n_inlet, n_outlet
    )

    GD = dict(
        Xno=Xno, conec=conec,
        n_inlet=n_inlet, n_outlet=n_outlet,
        T_solid_flat=T_solid_flat, k_e=k_e, k_n=k_n,
        T_arestas=T_arestas, mu_arestas=mu_arestas, C=C,
        A_net=A_net, A_net_sparse=A_net_sparse, np_net=np_net,
        K_mem=K_mem, M_mem=M_mem, h_ad=h_ad, nm=nm,
        Aglob=Aglob, U_sparse=U_sparse, pref=pref, vref=vref,
        params=params,
    )
    return GD


# ========================================================================
# 5. SIMULACAO TRANSIENTE DE VERIFICACAO DA BASE
# ========================================================================

def RodaSimulacaoBase(GD, t_max, p_inlet_dim):
    """
    Executa um laco temporal (Euler implicito) usando a matriz global
    do Gemeo Digital, com pressao de entrada constante `p_inlet_dim`
    [Pa], a fim de verificar a consistencia da base montada antes de
    iniciar os exercicios das Secoes 6.3 e 6.4.

    Retorna um dicionario `historico` com as series temporais de
    pressao no Outlet, vazao adimensional no Outlet, deslocamento no
    centro da membrana, alem dos campos finais (w, v) de toda a malha.
    """
    params = GD['params']
    dt = params['dt']
    nm = GD['nm']
    np_net = GD['np_net']
    n_inlet = GD['n_inlet']
    n_outlet = GD['n_outlet']
    pref = GD['pref']

    Nx_m, Ny_m = params['Nx_m'], params['Ny_m']
    no_centro = int((Ny_m // 2) * Nx_m + (Nx_m // 2))
    U_out = GD['U_sparse'][n_outlet, :].toarray().ravel()

    w_n = np.zeros(nm)
    v_n = np.zeros(nm)

    p_inlet_adim = p_inlet_dim / pref
    tempos = np.arange(0.0, t_max + 0.5 * dt, dt)

    hist_t, hist_p_out, hist_q_out, hist_w_centro = [], [], [], []

    for t in tempos:
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            GD['Aglob'], GD['M_mem'], w_n, v_n, p_inlet_adim, dt, nm, np_net, n_inlet
        )

        p_out_fisico = p_n[n_outlet] * pref
        q_out_adim = GD['h_ad'] ** 2 * np.sum(U_out * v_n)
        w_centro_fisico = w_n[no_centro] * (0.01 * params['R_fisico'])

        hist_t.append(t)
        hist_p_out.append(p_out_fisico)
        hist_q_out.append(q_out_adim)
        hist_w_centro.append(w_centro_fisico)

    historico = dict(
        t=np.array(hist_t), p_out=np.array(hist_p_out),
        q_out=np.array(hist_q_out), w_centro=np.array(hist_w_centro),
        w_final=w_n, v_final=v_n,
    )
    return historico


# ========================================================================
# 6. VISUALIZACAO DE DIAGNOSTICO DA BASE MONTADA
# ========================================================================

def PlotaEstadoBaseGD(GD, historico, save_path=None):
    """
    Gera um painel de diagnostico com tres paineis:
        (a) Mapa de temperatura na placa termica com a rede sobreposta;
        (b) Historico temporal de pressao e vazao no Outlet;
        (c) Forma final (deslocamento) da membrana elastica.

    Util apenas para validar visualmente a base montada antes de
    seguir para os exercicios das Secoes 6.3 e 6.4.
    """
    params = GD['params']
    Nx_t, Ny_t = params['Nx_t'], params['Ny_t']
    Lx, Ly = params['Lx'], params['Ly']

    fig = plt.figure(figsize=(16, 5))

    # --- (a) Campo de temperatura + rede sobreposta ---
    ax1 = fig.add_subplot(1, 3, 1)
    T_2D = GD['T_solid_flat'].reshape((Ny_t, Nx_t))
    x = np.linspace(0, Lx, Nx_t)
    y = np.linspace(0, Ly, Ny_t)
    cf = ax1.contourf(x, y, T_2D, levels=30, cmap='jet')
    for n1, n2 in GD['conec']:
        ax1.plot([GD['Xno'][n1, 0], GD['Xno'][n2, 0]],
                 [GD['Xno'][n1, 1], GD['Xno'][n2, 1]],
                 color='black', linewidth=0.4, alpha=0.6, zorder=2)
    ax1.set_aspect('equal')
    ax1.set_title('Placa Termica + Rede Hidraulica')
    ax1.set_xlabel('x (m)')
    ax1.set_ylabel('y (m)')
    fig.colorbar(cf, ax=ax1, label='T (°C)', shrink=0.8)

    # --- (b) Historico temporal de verificacao ---
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(historico['t'], historico['p_out'], 'r-', linewidth=1.5, label='Pressão Outlet (Pa)')
    ax2.set_xlabel('Tempo adimensional')
    ax2.set_ylabel('Pressão (Pa)', color='r')
    ax2.tick_params(axis='y', labelcolor='r')
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2b = ax2.twinx()
    ax2b.plot(historico['t'], historico['q_out'], 'b--', linewidth=1.2, label='Vazão Outlet (adim.)')
    ax2b.set_ylabel('Vazão adimensional', color='b')
    ax2b.tick_params(axis='y', labelcolor='b')
    ax2.set_title('Verificação Transiente da Base do GD')

    # --- (c) Forma final da membrana ---
    ax3 = fig.add_subplot(1, 3, 3, projection='3d')
    Nx_m, Ny_m = params['Nx_m'], params['Ny_m']
    Lx_ad, Ly_ad = params['Lx_ad'], params['Ly_ad']
    xm = np.linspace(0, Lx_ad, Nx_m)
    ym = np.linspace(0, Ly_ad, Ny_m)
    Xm, Ym = np.meshgrid(xm, ym)
    Wm = historico['w_final'].reshape((Ny_m, Nx_m))
    ax3.plot_surface(Xm, Ym, Wm, cmap='viridis', edgecolor='none', alpha=0.95)
    ax3.set_title('Forma Final da Membrana (adim.)')
    ax3.set_xlabel(r'$\hat{x}$')
    ax3.set_ylabel(r'$\hat{y}$')
    ax3.set_zlabel(r'$\hat{w}$')

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.show()
    return fig

# =======================================================================================
# EXERCÍCIOS - CAPÍTULO 6
# =======================================================================================

# ==============================================================================
# SECAO 6.3.2 -  Investigando o comportamento do sistema via Monte Carlo
# ==============================================================================

# ========================================================================
# SECAO 6.3.2 - EXERCICIO 1
# Análise Estacionária de Falhas Hidráulicas
# ========================================================================

def RandomFail(C_original, pO, fObs, rng):
    """
    Gera, de forma estocastica, uma instancia da rede com microcanais
    obstruidos.

    Cada canal possui probabilidade independente `pO` de estar
    obstruido. Quando um canal e obstruido, sua condutancia hidraulica
    nominal C_k e reduzida por um fator de severidade `fObs`, de modo
    que C_k <- C_k / fObs (conforme descrito na Secao 6.3 do PDF).

    Parameters
    ----------
    C_original : np.ndarray
        Vetor de condutancias nominais (sem falhas) de cada canal.
    pO : float
        Probabilidade individual de obstrucao de cada canal (0 a 1).
    fObs : float
        Fator de severidade da obstrucao (C_k e dividido por fObs).
    rng : np.random.Generator
        Gerador de numeros aleatorios (np.random.default_rng(...)).

    Returns
    -------
    C_modificado : np.ndarray
        Vetor de condutancias da realizacao estocastica gerada.
    """
    C_modificado = C_original.copy()
    obstruido = rng.random(C_original.shape[0]) < pO
    C_modificado[obstruido] = C_modificado[obstruido] / fObs
    return C_modificado


def ResolveRedeIsolada(conec, C, n_inlet, n_outlet, p_inlet, p_outlet):
    """
    Resolve o subsistema isolado (Rede Hidraulica + Placa Termica, ja
    embutida nas condutancias C) em regime permanente, com pressao
    prescrita em ambas as extremidades (Inlet e Outlet).

    A vazao total de entrada e calculada por "reacao", a partir da
    matriz de condutancia NAO restrita A (Equacao da Secao 6.3):

        q_inlet = sum_j A[n_inlet, j] * p_j   =>   qin = A[n_inlet,:] @ p
    """
    A = Assembly(conec, C)
    n = A.shape[0]
    Atilde = A.copy()
    b = np.zeros(n)

    Atilde[n_inlet, :] = 0.0
    Atilde[n_inlet, n_inlet] = 1.0
    b[n_inlet] = p_inlet

    Atilde[n_outlet, :] = 0.0
    Atilde[n_outlet, n_outlet] = 1.0
    b[n_outlet] = p_outlet

    pressao = np.linalg.solve(Atilde, b)
    q_inlet = A[n_inlet, :] @ pressao
    return q_inlet, pressao


def MonteCarloFalhaHidraulica(conec, C_nominal, n_inlet, n_outlet, p_inlet,
                               p_outlet, q_critico, pO, fObs, N, rng,
                               retornar_evolucao=False):
    """
    Executa N realizacoes de Monte Carlo do cenario de obstrucao
    estocastica dos microcanais (RandomFail), avaliando, em cada
    realizacao, se a vazao de entrada resultante fica abaixo do limite
    critico `q_critico`.

    Se `retornar_evolucao=True`, retorna o vetor da estimativa
    progressiva de Prob em funcao do numero de realizacoes (util para
    a analise de convergencia). Caso contrario, retorna apenas o valor
    final (escalar) da probabilidade estimada apos as N realizacoes.
    """
    eventos = np.empty(N, dtype=bool)
    for i in range(N):
        C_real = RandomFail(C_nominal, pO, fObs, rng)
        q_in, _ = ResolveRedeIsolada(conec, C_real, n_inlet, n_outlet, p_inlet, p_outlet)
        eventos[i] = q_in < q_critico

    if retornar_evolucao:
        return np.cumsum(eventos) / np.arange(1, N + 1)
    return float(np.mean(eventos))


def DeterminaTamanhoAmostralEstabilizacao(evolucao, margem_alvo=0.01, z_score=1.959964):
    """
    Define o tamanho amostral minimo N_min a partir do qual a estimativa
    progressiva de Monte Carlo `evolucao` (vetor Prob(N), N = 1, 2, ...)
    pode ser considerada estatisticamente estabilizada.

    CRITERIO TEORICO (Wald, aproximacao normal para uma proporcao): a
    estimativa Prob(N) e uma proporcao amostral (fracao de "sucessos"
    entre N realizacoes de Bernoulli), cujo erro padrao classico e

        SE(N) = sqrt( Phat * (1 - Phat) / N ),  Phat = evolucao[-1]

    com meia-largura do IC 95% dada por margem(N) = z_score * SE(N).
    N_min e o MENOR N a partir do qual margem(N) permanece, ate o final
    da amostragem, abaixo de `margem_alvo`.

    NOTA IMPORTANTE sobre o que esse criterio significa (e o que NAO
    significa): N_min e uma propriedade TEORICA de N e do valor final
    Phat - nao da trajetoria especifica realizada. Isso o torna
    deterministico e reprodutivel (nao depende da sorte da semente
    aleatoria), e e o que reproduz a ordem de grandeza de "algumas
    milhares de realizacoes" mencionada no PDF. Em contrapartida, ele
    NAO garante que a curva real, apos N_min, nunca mais oscile para
    fora da banda de confianca - isso e, alias, esperado: por
    definicao, um IC de 95% implica que, em ~5% dos pontos, a
    trajetoria realizada pode estar fora da banda mesmo apos a
    "estabilizacao" teorica. (Uma alternativa seria um criterio
    EMPIRICO, baseado em quando a trajetoria realizada de fato para de
    se afastar do valor final - mas esse criterio e estatisticamente
    fragil: depende de uma unica realizacao de Monte Carlo e pode
    variar muito de uma semente para outra, inclusive subestimando
    N_min por "sorte" da trajetoria, como observado para f_obs=10.)

    Parameters
    ----------
    evolucao : np.ndarray
        Estimativa progressiva de Prob, indices 0..N_total-1
        correspondendo a N = 1..N_total realizacoes.
    margem_alvo : float
        Meia-largura alvo do intervalo de confianca, em probabilidade
        (nao em %). Default: 0.01 (1 ponto percentual) - e o valor que
        reproduz a ordem de grandeza esperada pelo PDF. Requer
        N_convergencia grande o suficiente para que a margem de 1 p.p.
        seja de fato atingida dentro da amostra simulada (caso
        contrario, N_min fica preso no limite superior da amostra).
    z_score : float
        Quantil normal associado ao nivel de confianca desejado.
        Default: 1.959964 (95% de confianca).

    Returns
    -------
    N_min : int
        Tamanho amostral minimo para estabilizacao, segundo o criterio
        teorico acima.
    margem_final : float
        Margem teorica do IC 95% em N_min (deve ser <= margem_alvo, e
        permanece abaixo disso ate o final da amostragem).
    margem_vs_N : np.ndarray
        Vetor completo margem(N) para N = 1..N_total, usado tambem para
        plotar a banda de confianca (funil) em torno do valor
        convergido.
    """
    N_total = len(evolucao)
    Ns = np.arange(1, N_total + 1)
    phat_final = evolucao[-1]
    margem_vs_N = z_score * np.sqrt(phat_final * (1.0 - phat_final) / Ns)

    dentro_da_margem = margem_vs_N <= margem_alvo

    idx_estab = N_total - 1
    for i in range(N_total - 1, -1, -1):
        if dentro_da_margem[i]:
            idx_estab = i
        else:
            break

    N_min = idx_estab + 1   # indice 0 corresponde a N = 1 realizacao
    margem_final = margem_vs_N[idx_estab]
    return N_min, margem_final, margem_vs_N


def RodaExercicio1_FalhasHidraulicas(GD, q_critico=1.25e-5,
                                      N_convergencia=5000,
                                      pO_convergencia=0.6,
                                      fObs_lista=(5, 10),
                                      pO_grid=None,
                                      N_final=3000,
                                      p_outlet=0.0,
                                      semente=42):
    """
    Resolve o Exercicio 1 da Secao 6.3.2 do PDF (Analise Estacionaria
    de Falhas Hidraulicas):

        (a) Implementa RandomFail() para gerar obstrucoes estocasticas
            das condutancias dos canais;
        (b) Estuda a convergencia da estimativa de Monte Carlo Prob(N)
            para um cenario de obstrucao representativo
            (p_O = `pO_convergencia`), permitindo definir um tamanho
            amostral N minimo para a estabilizacao estatistica;
        (c) Varre o dominio p_O em [0.05, 0.65] e calcula os valores
            CONVERGIDOS de Prob (usando N = `N_final` amostras por
            ponto), contrastando dois fatores de severidade de
            obstrucao distintos: f_obs = 5 e f_obs = 10.

    A rede e a placa termica utilizadas sao exatamente as do Gemeo
    Digital nominal montado em `GD` (MontaGemeoDigital): as
    condutancias `GD['C']` ja incorporam o efeito da temperatura local
    de cada canal (subsistema isolado Rede + Placa Termica).

    Retorna um dicionario com todos os resultados numericos gerados.
    """
    if pO_grid is None:
        pO_grid = np.arange(0.05, 0.65 + 1e-9, 0.05)

    conec = GD['conec']
    C_nominal = GD['C']
    n_inlet = GD['n_inlet']
    n_outlet = GD['n_outlet']
    p_inlet = GD['params']['p_inlet_dim']   # 5000 Pa nominal (Secao 6.2)

    # ---- (0) Diagnostico de calibracao: vazao NOMINAL (sem falhas) ----
    # Compara a vazao de entrada do sistema isolado, sem nenhuma
    # obstrucao, com o limite critico (usado apenas internamente, para
    # compor `resultados['razao_nominal']`; impressao no terminal
    # removida a pedido).
    q_nominal, _ = ResolveRedeIsolada(conec, C_nominal, n_inlet, n_outlet,
                                       p_inlet, p_outlet)
    razao_nominal = q_nominal / q_critico

    # ---- (b) Estudo de convergencia da estimativa de Monte Carlo ------
    print(f"\n[Exercicio 1] Estudo de convergencia em p_O = {pO_convergencia} "
          f"para N = {N_convergencia} realizacoes...")
    evolucoes_convergencia = {}
    N_estabilizacao = {}
    margem_estabilizacao = {}
    margem_vs_N_convergencia = {}
    for fObs in fObs_lista:
        rng_conv = np.random.default_rng(semente)
        evol = MonteCarloFalhaHidraulica(
            conec, C_nominal, n_inlet, n_outlet, p_inlet, p_outlet,
            q_critico, pO_convergencia, fObs, N_convergencia,
            rng=rng_conv, retornar_evolucao=True
        )
        evolucoes_convergencia[fObs] = evol

        N_min, margem, margem_vs_N = DeterminaTamanhoAmostralEstabilizacao(evol)
        N_estabilizacao[fObs] = N_min
        margem_estabilizacao[fObs] = margem
        margem_vs_N_convergencia[fObs] = margem_vs_N

        print(f"   f_obs={fObs:<3d} -> Prob estabilizada (N={N_convergencia}): "
              f"{evol[-1]*100:.2f}%")
        print(f"             -> Tamanho amostral minimo para estabilizacao "
              f"estatistica: N_min = {N_min} realizacoes "
              f"(margem teorica do IC 95% da proporcao: ±{margem*100:.2f} p.p.)")
        if N_min >= N_convergencia:
            print(f"             -> AVISO: a margem alvo so foi atingida no "
                  f"limite da amostra simulada (N={N_convergencia}); "
                  f"considere aumentar N_convergencia para confirmar a "
                  f"estabilizacao plena.")

    N_simulacoes_convergencia = N_convergencia * len(fObs_lista)

    # ---- (c) Varredura em p_O para os dois fatores de severidade ------
    print(f"\n[Exercicio 1] Varredura em p_O em [{pO_grid[0]:.2f}, {pO_grid[-1]:.2f}] "
          f"com N = {N_final} realizacoes por ponto...")
    prob_vs_pO = {fObs: np.zeros_like(pO_grid) for fObs in fObs_lista}
    rng_varredura = np.random.default_rng(semente + 1)

    for fObs in fObs_lista:
        for idx, pO in enumerate(pO_grid):
            prob = MonteCarloFalhaHidraulica(
                conec, C_nominal, n_inlet, n_outlet, p_inlet, p_outlet,
                q_critico, pO, fObs, N_final,
                rng=rng_varredura, retornar_evolucao=False
            )
            prob_vs_pO[fObs][idx] = prob
        print(f"   f_obs={fObs:<3d} concluido.")

    N_simulacoes_varredura = N_final * len(pO_grid) * len(fObs_lista)

    N_simulacoes_total = N_simulacoes_convergencia + N_simulacoes_varredura

    resultados = dict(
        q_critico=q_critico,
        q_nominal=q_nominal,
        razao_nominal=razao_nominal,
        p_inlet=p_inlet, p_outlet=p_outlet,
        N_convergencia=N_convergencia,
        pO_convergencia=pO_convergencia,
        evolucoes_convergencia=evolucoes_convergencia,
        N_estabilizacao=N_estabilizacao,
        margem_estabilizacao=margem_estabilizacao,
        margem_vs_N_convergencia=margem_vs_N_convergencia,
        pO_grid=pO_grid,
        N_final=N_final,
        prob_vs_pO=prob_vs_pO,
        fObs_lista=fObs_lista,
        N_simulacoes_convergencia=N_simulacoes_convergencia,
        N_simulacoes_varredura=N_simulacoes_varredura,
        N_simulacoes_total=N_simulacoes_total,
    )
    return resultados


def PlotaExercicio1_FalhasHidraulicas(resultados,
                                       save_path_convergencia=None,
                                       save_path_varredura=None):
    """
    Gera os dois graficos pedidos no Exercicio 1 da Secao 6.3.2:

        (1) Evolucao da estimativa de Prob em funcao do numero
            progressivo de realizacoes de Monte Carlo, permitindo
            avaliar o comportamento assintotico/estabilizacao;
        (2) Valores CONVERGIDOS de Prob em funcao da probabilidade de
            obstrucao individual p_O, contrastando f_obs = 5 e
            f_obs = 10.
    """
    # --- Grafico 1: convergencia da estimativa Prob(N) -----------------
    fig1, ax1 = plt.subplots(figsize=(9, 5.5))
    N_conv = resultados['N_convergencia']
    eixo_N = np.arange(1, N_conv + 1)
    cores = plt.rcParams['axes.prop_cycle'].by_key()['color']
    N_estab_dict = resultados.get('N_estabilizacao', {})
    margem_dict = resultados.get('margem_estabilizacao', {})
    margem_vs_N_dict = resultados.get('margem_vs_N_convergencia', {})

    for i, (fObs, evol) in enumerate(resultados['evolucoes_convergencia'].items()):
        cor = cores[i % len(cores)]
        valor_final = evol[-1]

        # --- Banda de confianca em "funil" (IC 95%, estreita com N) --------
        # IMPORTANTE: a banda e centrada no VALOR FINAL CONVERGIDO (constante),
        # nao na curva ruidosa `evol`. Centralizar em `evol` faz a banda
        # "tremer" junto com as flutuacoes da trajetoria de Monte Carlo,
        # destruindo o formato de funil liso esperado teoricamente
        # (a banda deve refletir apenas a evolucao assintotica de 1/sqrt(N)
        # em torno do valor para o qual o estimador converge).
        if fObs in margem_vs_N_dict:
            limite_inf = np.clip(valor_final - margem_vs_N_dict[fObs], 0.0, 1.0)
            limite_sup = np.clip(valor_final + margem_vs_N_dict[fObs], 0.0, 1.0)
            ax1.fill_between(eixo_N, limite_inf, limite_sup,
                              color=cor, alpha=0.15, linewidth=0,
                              label=fr'Banda IC 95% ($f_{{obs}}$={fObs})')

        # --- Curva da estimativa progressiva --------------------------------
        ax1.plot(eixo_N, evol, linewidth=1.3, color=cor,
                 label=fr'Estimativa de Prob ($f_{{obs}}$ = {fObs})')

        # --- Linha pontilhada no valor convergido ---------------------------
        ax1.axhline(valor_final, color=cor, linestyle='--', linewidth=1.3,
                    alpha=0.9,
                    label=fr'Valor convergido $f_{{obs}}$={fObs} $\approx$ {valor_final:.4f}')

        # --- Marca o N minimo de estabilizacao ------------------------------
        if fObs in N_estab_dict:
            N_min = N_estab_dict[fObs]
            ax1.plot(N_min, evol[N_min - 1], 'o', color=cor,
                      markersize=7, markeredgecolor='black', markeredgewidth=0.8,
                      zorder=5)
            ax1.annotate(fr'$N_{{min}}$={N_min}',
                         xy=(N_min, evol[N_min - 1]),
                         xytext=(8, 10 if i % 2 == 0 else -18),
                         textcoords='offset points', fontsize=9, color=cor,
                         fontweight='bold')

    ax1.set_xscale('log')
    ax1.set_xlim(left=10, right=N_conv)
    ax1.set_ylim(0.8, 1)
    ax1.set_xlabel('Número de realizações de Monte Carlo, N')
    ax1.set_ylabel(r'Estimativa progressiva de $Prob(q_{inlet} < q_{crit})$')
    ax1.set_title(fr"Convergência da estimativa de Monte Carlo ($p_O$ = {resultados['pO_convergencia']})"
                  "\n(banda = IC 95% teórico da proporção; linha tracejada = valor convergido; "
                  "marcador = $N_{min}$ teórico)")
    ax1.grid(True, which='both', linestyle=':', alpha=0.6)
    ax1.legend(fontsize=8, ncol=1)
    plt.tight_layout()
    if save_path_convergencia:
        plt.savefig(save_path_convergencia, dpi=200, bbox_inches='tight')
    plt.show()

    # --- Grafico 2: Prob convergida vs p_O ------------------------------
    fig2, ax2 = plt.subplots(figsize=(8, 5.5))
    estilos = {5: 'o-', 10: 's--'}
    for fObs in resultados['fObs_lista']:
        ax2.plot(resultados['pO_grid'], resultados['prob_vs_pO'][fObs] * 100.0,
                 estilos.get(fObs, '.-'), linewidth=1.6, markersize=5,
                 label=fr'$f_{{obs}}$ = {fObs}')
    ax2.set_xlabel(r'Probabilidade de obstrução individual, $p_O$')
    ax2.set_ylabel(r'Probabilidade estimada, $Prob(q_{inlet} < q_{crit})$ (%)')
    ax2.set_title(fr"Probabilidade de Falha Hidráulica vs $p_O$  (N = {resultados['N_final']})")
    ax2.grid(True, linestyle=':', alpha=0.6)
    ax2.legend()
    plt.tight_layout()
    if save_path_varredura:
        plt.savefig(save_path_varredura, dpi=200, bbox_inches='tight')
    plt.show()

    return fig1, fig2

# ==============================================================================
# SEÇÃO 6.3.2 - EXERCÍCIO 2
# Análise Dinâmica do Gêmeo Digital Completo
# ==============================================================================

def MontaSistemaFalho(GD, C_real, dt):
    """
    Reconstroi, para uma instancia de rede COM FALHAS (vetor de
    condutancias `C_real`, ja gerado por RandomFail) e um passo de
    tempo `dt`, a matriz global acoplada Aglob (necessaria para marchar
    no tempo via Resolve_Passo_Tempo) E a matriz de condutancia
    hidraulica "fisica" 𝔸 = D^T K D, adimensionalizada e SEM a
    substituicao da linha do Inlet, usada na forma quadratica da
    potencia instantanea P(t) = p(t)^T 𝔸 p(t).

    Nota: como o vetor de condutancias muda a cada realizacao
    estocastica (e o passo de tempo pode mudar entre cenarios, ja que o
    Exercicio 2 pede dt=0.05 e dt=0.1), tanto a rede (A_net) quanto
    Aglob (que depende de dt) precisam ser remontados aqui -- somente
    as matrizes da membrana (K_mem, M_mem) sao reaproveitadas de GD,
    pois nao dependem da falha hidraulica.

    Retorna (Aglob_falha, A_potencia, pref).
    """
    params = GD['params']
    conec = GD['conec']

    # ---- Rede hidraulica com falhas, para esta realizacao --------------
    A_net_falha = Assembly(conec, C_real)
    A_net_falha_sparse = sparse.csr_matrix(A_net_falha)

    # ---- Matriz global acoplada (membrana + rede), para este dt -------
    Aglob_falha, _U_sparse, pref, vref = Monta_Matriz_Global_Acoplada(
        GD['K_mem'], GD['M_mem'], A_net_falha_sparse,
        params['Nx_m'], params['Ny_m'], params['Lx_ad'], params['Ly_ad'],
        params['R_ad'], GD['h_ad'], dt, params['beta_ad'],
        params['sigma'], params['rho'], params['e_esp'], params['R_fisico'],
        GD['n_inlet'], GD['n_outlet']
    )

    # ---- Matriz "fisica" D^T K D, adimensionalizada, SEM a -------------
    # substituicao da linha do Inlet (usada apenas para o calculo de
    # P(t), e nao para a resolucao do sistema). O fator de escala e o
    # mesmo aplicado internamente na montagem de Aglob (Secao 5.2.4 do
    # PDF: A_adim = A_dimensional * pref / (vref * R_fisico**2)).
    fator_escala = pref / (vref * params['R_fisico'] ** 2)
    A_potencia = A_net_falha_sparse * fator_escala

    return Aglob_falha, A_potencia, pref


def SimulaEnergiaDissipada(GD, C_real, dt, t_max, p_inlet_dim, t_i=0.0):
    """
    Executa UMA realizacao transiente completa do Gemeo Digital
    (acoplamento forte multifisico: membrana + rede + placa, ja
    embutida em GD['C']/`C_real`), para uma instancia de rede com
    falhas `C_real`, e calcula o indicador energetico

        E = integral_{t_i}^{t_f} P(t) dt ,   P(t) = p(t)^T 𝔸 p(t)

    via a regra do trapezio (np.trapz), integrando a potencia
    hidraulica instantanea dissipada no reservatorio ao longo da
    trajetoria temporal.

    Retorna (E, P_hist, tempos).
    """
    nm = GD['nm']
    np_net = GD['np_net']
    n_inlet = GD['n_inlet']

    Aglob_falha, A_potencia, pref = MontaSistemaFalho(GD, C_real, dt)

    p_inlet_adim = p_inlet_dim / pref

    w_n = np.zeros(nm)
    v_n = np.zeros(nm)

    tempos = np.arange(t_i, t_max + 0.5 * dt, dt)
    P_hist = np.empty(tempos.shape[0])

    for i in range(tempos.shape[0]):
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob_falha, GD['M_mem'], w_n, v_n, p_inlet_adim, dt, nm, np_net, n_inlet
        )
        P_hist[i] = p_n @ (A_potencia @ p_n)

    E = _integral_trapezoidal(P_hist, dx=dt)
    return E, P_hist, tempos


def MonteCarloEnergiaDinamica(GD, pO, fObs, dt, N, t_max=4.0,
                               p_inlet_dim=None, E_critico=7.0,
                               rng=None, retornar_amostras=False,
                               imprimir_progresso=True):
    """
    Executa N realizacoes de Monte Carlo do GD COMPLETO (acoplamento
    forte multifisico, Secao 6.3.2, Item 2), cada uma com uma instancia
    de falha estocastica distinta dos microcanais (RandomFail, mesma
    hipotese p_O / f_obs do Item 1), avaliando se a energia total
    dissipada E fica abaixo do limite critico `E_critico`.

    Retorna a probabilidade estimada Prob(E < E_critico) e,
    opcionalmente (retornar_amostras=True), o vetor completo das N
    amostras de E (util para histogramas e analises adicionais).
    """
    if p_inlet_dim is None:
        p_inlet_dim = GD['params']['p_inlet_dim']
    if rng is None:
        rng = np.random.default_rng()

    C_nominal = GD['C']
    E_amostras = np.empty(N)

    t0 = time.time()
    passo_print = max(1, N // 20)
    for i in range(N):
        C_real = RandomFail(C_nominal, pO, fObs, rng)
        E_amostras[i], _, _ = SimulaEnergiaDissipada(GD, C_real, dt, t_max, p_inlet_dim)
        if imprimir_progresso and (i + 1) % passo_print == 0:
            decorrido = time.time() - t0
            print(f"      ... {i + 1}/{N} realizacoes (f_obs={fObs}, dt={dt}) "
                  f"- {decorrido:.1f}s decorridos")

    prob = float(np.mean(E_amostras < E_critico))
    if retornar_amostras:
        return prob, E_amostras
    return prob


def RodaExercicio2_AnaliseDinamica(GD, E_critico=7.0, pO=0.6,
                                    fObs_lista=(5, 10), dt_lista=(0.05, 0.1),
                                    N=2000, t_max=4.0, p_inlet_dim=None,
                                    semente=123):
    """
    Resolve o Exercicio 2 da Secao 6.3.2 do PDF (Analise Dinamica do
    Gemeo Digital Completo):

        - Acoplamento forte multifisico completo (membrana 51x51 +
          rede hidraulica + placa termica, ja embutidos em `GD`);
        - Falha estocastica dos microcanais sob probabilidade
          individual p_O (mesma RandomFail() do Item 1);
        - Indicador energetico E = integral_0^4 P(t) dt, com
          P(t) = p(t)^T 𝔸 p(t);
        - Estima Prob(E < E_critico) para dt = 0.05 e dt = 0.1, com
          N = 2000 realizacoes cada, contrastando dois fatores de
          severidade de obstrucao f_obs = 5 e f_obs = 10 (mesmos
          valores investigados no Exercicio 1).

    NOTA SOBRE p_O: o enunciado do Item 2 nao especifica um novo valor
    para p_O (apenas diz "mantendo a hipotese de falha sob
    probabilidade individual p_O"), entao adota-se aqui, por padrao, o
    MESMO cenario representativo ja usado no estudo de convergencia do
    Exercicio 1 (p_O = 0.6). Ajuste o argumento `pO` se um valor
    diferente for desejado.

    Retorna um dicionario com todos os resultados numericos gerados,
    indexados por (f_obs, dt).
    """
    if p_inlet_dim is None:
        p_inlet_dim = GD['params']['p_inlet_dim']

    prob_dict = {}
    amostras_dict = {}

    contador_semente = 0
    for fObs in fObs_lista:
        for dt in dt_lista:
            print(f"\n[Exercicio 2] Rodando MC dinamico: f_obs={fObs}, dt={dt}, "
                  f"N={N} realizacoes (p_O={pO}) ...")
            rng = np.random.default_rng(semente + contador_semente)
            contador_semente += 1

            prob, E_amostras = MonteCarloEnergiaDinamica(
                GD, pO, fObs, dt, N, t_max=t_max, p_inlet_dim=p_inlet_dim,
                E_critico=E_critico, rng=rng, retornar_amostras=True
            )
            prob_dict[(fObs, dt)] = prob
            amostras_dict[(fObs, dt)] = E_amostras
            print(f"   -> Prob(E < {E_critico}) = {prob * 100:.2f}%   "
                  f"(E_medio={E_amostras.mean():.3f}, "
                  f"E_min={E_amostras.min():.3f}, E_max={E_amostras.max():.3f})")

    resultados = dict(
        E_critico=E_critico, pO=pO, fObs_lista=fObs_lista,
        dt_lista=dt_lista, N=N, t_max=t_max,
        prob=prob_dict, amostras=amostras_dict,
    )
    return resultados

def PlotaExercicio2_AnaliseDinamica(resultados, save_path=None):
    """
    Gera um painel com dois graficos comparativos para o Exercicio 2:

        (a) Prob(E < E_critico) em funcao do passo de tempo dt,
            contrastando os dois fatores de severidade f_obs;
        (b) Histogramas sobrepostos da distribuicao amostral de E para
            cada combinacao (f_obs, dt), permitindo visualizar o
            impacto do passo de tempo na precisao estatistica do
            estimador (dispersao da amostra em torno do limite
            E_critico).
    """
    fObs_lista = resultados['fObs_lista']
    dt_lista = resultados['dt_lista']
    E_critico = resultados['E_critico']

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    cores = plt.rcParams['axes.prop_cycle'].by_key()['color']

    # --- (a) Prob(E<E_critico) vs dt, para cada f_obs --------------------
    ax1 = axes[0]
    for i, fObs in enumerate(fObs_lista):
        probs = [resultados['prob'][(fObs, dt)] * 100.0 for dt in dt_lista]
        ax1.plot(dt_lista, probs, 'o-', color=cores[i % len(cores)],
                  linewidth=1.8, markersize=8, label=fr'$f_{{obs}}$ = {fObs}')
    ax1.set_xlabel(r'Passo de tempo, $\delta t$')
    ax1.set_ylabel(fr'Prob($E < {E_critico}$)  (%)')
    ax1.set_title("Análise Dinâmica do GD Completo\n"
                  fr"($p_O$={resultados['pO']}, N={resultados['N']} por ponto)")
    ax1.set_xticks(dt_lista)
    ax1.grid(True, linestyle=':', alpha=0.6)
    ax1.legend()

    # --- (b) Histograma da energia dissipada E, por combinacao ----------
    ax2 = axes[1]
    estilos_dt = ['-', '--', '-.', ':']
    for i, fObs in enumerate(fObs_lista):
        for j, dt in enumerate(dt_lista):
            E_amostras = resultados['amostras'][(fObs, dt)]
            ax2.hist(E_amostras, bins=40, histtype='step', linewidth=1.6,
                      color=cores[i % len(cores)],
                      linestyle=estilos_dt[j % len(estilos_dt)],
                      label=fr'$f_{{obs}}$={fObs}, $\delta t$={dt}')
    ax2.axvline(E_critico, color='black', linestyle=':', linewidth=1.8,
                label=fr'$E_{{crítico}}$ = {E_critico}')
    ax2.set_xlabel('Energia dissipada total, E')
    ax2.set_ylabel('Frequência (contagem de realizações)')
    ax2.set_title('Distribuição amostral de E')
    ax2.grid(True, linestyle=':', alpha=0.5)
    ax2.legend(fontsize=8)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=200, bbox_inches='tight')
    plt.show()
    return fig


# ==============================================================================
# SECAO 6.4.3 Investigando o comportamento do sistema via aproximação de dados
# ==============================================================================

# ========================================================================
# SEÇÃO 6.4.3 - EXERCÍCIO 1
# Interpolação Linear Local
# ========================================================================

def Interpola_Potencia_Linear(t_dados, P_dados):
    """ 
    Exercício 1: Interpolação Linear Local (Splines lineares).
    Retorna uma função avaliável que liga os dados por retas.
    """
    return interp1d(t_dados, P_dados, kind='linear')


# ========================================================================
# SEÇÃO 6.4.3 - EXERCÍCIO 2
# Interpolação Cúbica Local
# ========================================================================

def Interpola_Potencia_Cubica(t_dados, P_dados):
    """ 
    Exercício 2: Interpolação Cúbica Local (Splines cúbicos).
    Retorna uma função avaliável garantindo suavidade nas interfaces.
    """
    return interp1d(t_dados, P_dados, kind='cubic')


# ========================================================================
# # SEÇÃO 6.4.3 - EXERCÍCIO 3
# Regressão Polinomial Global
# ========================================================================



# ========================================================================
# # SEÇÃO 6.4.3 - EXERCÍCIO 3
# SEÇÃO 6.4.3 - EXERCÍCIO 3
# Regressão Polinomial Global
# ========================================================================

# parte da nat vai aqui
def Regressao_Polinomial_Global(t_dados, P_dados, graus=None):
    """
    Exercício 3: Regressão Polinomial Global pelo Método dos Mínimos Quadrados.

    Para cada grau m em `graus`, ajusta um polinômio global aos dados
    (t_dados, P_dados) e calcula o erro L2 (Eq. 6.14):
        e = sqrt( integral_{t0}^{tf} [p(t) - P(t)]^2 dt )
    via regra do trapézio.

    Retorna um dicionário com:
        'graus'        : lista de graus testados
        'coeficientes' : {m: coefs} — coeficientes do np.polyfit
        'erros_L2'     : {m: e}     — erro L2 de cada grau
    """
    if graus is None:
        graus = list(range(3, 16))

    coeficientes = {}
    erros_L2 = {}

    for m in graus:
        coefs = np.polyfit(t_dados, P_dados, m)
        p_aprox = np.polyval(coefs, t_dados)
        residuo_quad = (p_aprox - P_dados) ** 2
        erro = np.sqrt(_integral_trapezoidal(residuo_quad,
                                             x=t_dados))
        coeficientes[m] = coefs
        erros_L2[m] = erro

    return {'graus': graus, 'coeficientes': coeficientes, 'erros_L2': erros_L2}


def SimulaEnergiaDissipada_Ruidosa(GD, dt, t_max, p_inlet_dim, seed=None):
    """
    Exercício 4: variante ruidosa de SimulaEnergiaDissipada.

    Em cada passo de tempo aplica uma perturbação uniforme U(-0.15, 0.15)
    à pressão de entrada adimensional:
        p_inlet(t) = (p_inlet_dim / p_ref) * [1 + U(-0.15, 0.15)]

    Retorna (E, P_hist, tempos) — mesma assinatura que
    SimulaEnergiaDissipada, porém com sinal de potência ruidoso.
    """
    rng = np.random.default_rng(seed)

    nm = GD['nm']
    np_net = GD['np_net']
    n_inlet = GD['n_inlet']

    Aglob, A_potencia, pref = MontaSistemaFalho(GD, GD['C'], dt)

    p_inlet_adim_base = p_inlet_dim / pref

    w_n = np.zeros(nm)
    v_n = np.zeros(nm)

    tempos = np.arange(0.0, t_max + 0.5 * dt, dt)
    P_hist = np.empty(tempos.shape[0])

    for i in range(tempos.shape[0]):
        ruido = -0.15 + 0.3 * rng.random()
        p_inlet_adim = p_inlet_adim_base * (1.0 + ruido)
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob, GD['M_mem'], w_n, v_n, p_inlet_adim, dt, nm, np_net, n_inlet
        )
        P_hist[i] = p_n @ (A_potencia @ p_n)

    E = _integral_trapezoidal(P_hist, dx=dt)
    return E, P_hist, tempos


# ========================================================================
# SEÇÃO 6.4.3 - EXERCÍCIO 4
# Análise de Sensibilidade a Ruídos Estocásticos
# ========================================================================


def Analise_Ruido_Estocastico(GD, dt, t_max, p_inlet_dim,
                               graus=None, N_realizacoes=10, seed=42):
    """
    Exercício 4: Análise de sensibilidade a ruídos estocásticos.

    Gera N_realizacoes trajetórias ruidosas de P(t) e, para cada uma,
    aplica a regressão polinomial global (Exercício 3) nos graus
    especificados.

    Retorna um dicionário com:
        'tempos'          : vetor de tempo (igual para todas as realizações)
        'P_ruidoso_medio' : média de P(t) sobre as N realizações
        'P_realizacoes'   : lista com os N vetores P_hist ruidosos
        'erros_L2_medio'  : {m: média do erro L2 sobre as N realizações}
        'erros_L2_todas'  : {m: lista com N erros L2}
        'graus'           : graus testados
    """
    if graus is None:
        graus = list(range(3, 16))

    rng_mestre = np.random.default_rng(seed)

    P_realizacoes = []
    tempos = None

    for k in range(N_realizacoes):
        s = int(rng_mestre.integers(0, 2**31))
        _E, P_hist, t_vec = SimulaEnergiaDissipada_Ruidosa(
            GD, dt, t_max, p_inlet_dim, seed=s
        )
        P_realizacoes.append(P_hist)
        if tempos is None:
            tempos = t_vec

    P_stack = np.array(P_realizacoes)          # (N_realizacoes, n_steps)
    P_ruidoso_medio = P_stack.mean(axis=0)

    erros_L2_todas = {m: [] for m in graus}
    for k in range(N_realizacoes):
        res = Regressao_Polinomial_Global(tempos, P_realizacoes[k], graus=graus)
        for m in graus:
            erros_L2_todas[m].append(res['erros_L2'][m])

    erros_L2_medio = {m: float(np.mean(erros_L2_todas[m])) for m in graus}

    return {
        'tempos': tempos,
        'P_ruidoso_medio': P_ruidoso_medio,
        'P_realizacoes': P_realizacoes,
        'erros_L2_medio': erros_L2_medio,
        'erros_L2_todas': erros_L2_todas,
        'graus': graus,
    }

# ==============================================================================
# SEÇÃO 6.4.5 Investigando o comportamento do sistema via diferenciação numérica
# ==============================================================================

# ========================================================================
# SEÇÃO 6.4.5 - EXERCÍCIO 1
# Análise Numérica de Sensibilidade
# ========================================================================

def CopiaParams645(params_base, TC=None, H_um=None, dt=None):
    """
    Copia os parâmetros do caso nominal e altera apenas TC, H ou dt.

    TC   em graus Celsius.
    H_um em micrometros.
    """
    params = params_base.copy()

    if TC is not None:
        params["TC"] = float(TC)

    if H_um is not None:
        params["H_k"] = float(H_um) * 1e-6

    if dt is not None:
        params["dt"] = float(dt)

    return params


def _matriz_A_scaled_sem_bc(GD):
    """
    Reconstrói a matriz hidráulica adimensional SEM impor a linha de pressão
    prescrita no inlet.

    Essa matriz é usada para:
        P(t) = p^T A p
        q_inlet = linha_inlet(A) @ p
    """
    params = GD["params"]
    fator_escala = GD["pref"] / (GD["vref"] * params["R_fisico"] ** 2)
    return sparse.csr_matrix(GD["A_net"]) * fator_escala

def _row_dot_scalar(A, i, x):
    """
    Calcula produto da linha i da matriz A pelo vetor x e devolve escalar float.
    Funciona mesmo quando A é matriz esparsa e o resultado sai como array 1D/2D.
    """
    return float(np.asarray(A.getrow(i).dot(x)).ravel()[0])


def _derivada_A_scaled_H_sem_bc(GD):
    """
    Como A(H) = H^4 A0, então dA/dH = (4/H) A(H).

    Aqui H está em metros, então essa derivada é em relação a H_m.
    """
    H_m = GD["params"]["H_k"]
    A_scaled = _matriz_A_scaled_sem_bc(GD)
    return (4.0 / H_m) * A_scaled


def _derivada_A_scaled_H_com_bc(GD):
    """
    Derivada da matriz hidráulica adimensional COM a condição de pressão
    prescrita no inlet.

    A linha do inlet vira uma linha da identidade em Monta_Matriz_Global_Acoplada.
    Como essa linha não depende de H, sua derivada é zero.
    """
    dA = _derivada_A_scaled_H_sem_bc(GD).tolil()
    n_inlet = GD["n_inlet"]
    dA[n_inlet, :] = 0.0
    return dA.tocsr()


def SimulaCaso645(params_base, TC=None, H_um=None, t_max=4.0, dt=0.05,
                  calcular_direto_H=False):
    """
    Roda uma simulação transiente do Gêmeo Digital para um valor de TC e H.

    Retorna:
        E       = integral de P(t) no intervalo [0, t_max]
        qin_tf  = vazão de entrada no tempo final
        V_tf    = volume acumulado no reservatório no tempo final

    Se calcular_direto_H=True, também calcula as derivadas diretas:
        dE_dH_direct
        dqin_dH_direct
        dV_dH_direct

    As derivadas diretas são retornadas em relação a H em micrometros,
    para ficarem comparáveis com as diferenças finitas feitas em H_um.
    """
    params = CopiaParams645(params_base, TC=TC, H_um=H_um, dt=dt)
    GD = MontaGemeoDigital(params)

    nm = GD["nm"]
    np_net = GD["np_net"]
    n_inlet = GD["n_inlet"]
    pref = GD["pref"]
    vref = GD["vref"]
    h_ad = GD["h_ad"]

    p_inlet_dim = params["p_inlet_dim"]
    p_inlet_adim = p_inlet_dim / pref

    Aglob = GD["Aglob"]
    M_mem = GD["M_mem"]

    A_scaled = _matriz_A_scaled_sem_bc(GD)
    qref = params["R_fisico"] ** 2 * vref
    wref = 0.01 * params["R_fisico"]
    area_elemento_fisico = (h_ad * params["R_fisico"]) ** 2

    if calcular_direto_H:
        dA_scaled_sem_bc_dH_m = _derivada_A_scaled_H_sem_bc(GD)
        dA_scaled_com_bc_dH_m = _derivada_A_scaled_H_com_bc(GD)

    w_n = np.zeros(nm)
    v_n = np.zeros(nm)

    sw_n = np.zeros(nm)
    sv_n = np.zeros(nm)
    sp_n = np.zeros(np_net)

    tempos = np.arange(0.0, t_max + 0.5 * dt, dt)

    P_hist = []
    qin_hist = []
    V_hist = []

    dP_dH_hist_m = []

    for _t in tempos:
        # Estado principal: resolve [w, v, p]
        w_n, v_n, p_n = Resolve_Passo_Tempo(
            Aglob, M_mem, w_n, v_n, p_inlet_adim,
            dt, nm, np_net, n_inlet
        )

        # Potência adimensional/coerente com o exercício 6.3.2
        P_inst = float(p_n @ (A_scaled @ p_n))

        # Vazão de entrada por reação hidráulica
        qin_adim = _row_dot_scalar(A_scaled, n_inlet, p_n)
        qin_fis = qin_adim * qref

        # Volume físico acumulado no reservatório
        V_fis = float(np.sum(w_n * wref) * area_elemento_fisico)

        P_hist.append(P_inst)
        qin_hist.append(qin_fis)
        V_hist.append(V_fis)

        if calcular_direto_H:
            # Sistema das sensibilidades:
            # G [sw, sv, sp] = [sw_n/dt, M sv_n/dt, -dA/dH p]
            b_sw = sw_n / dt
            b_sv = M_mem.dot(sv_n) / dt
            b_sp = -(dA_scaled_com_bc_dH_m @ p_n)

            b_sens = np.concatenate([b_sw, b_sv, b_sp])
            sol_sens = spsolve(Aglob, b_sens)

            sw_n = sol_sens[0:nm]
            sv_n = sol_sens[nm:2 * nm]
            sp_n = sol_sens[2 * nm:]

            dP_dH_m = float(
                2.0 * (sp_n @ (A_scaled @ p_n))
                + p_n @ (dA_scaled_sem_bc_dH_m @ p_n)
            )
            dP_dH_hist_m.append(dP_dH_m)

    P_hist = np.array(P_hist)
    qin_hist = np.array(qin_hist)
    V_hist = np.array(V_hist)

    E = float(_integral_trapezoidal(P_hist, dx=dt))

    saida = {
        "TC": params["TC"],
        "H_um": params["H_k"] * 1e6,
        "t": tempos,
        "P": P_hist,
        "qin": qin_hist,
        "V": V_hist,
        "E": E,
        "qin_tf": float(qin_hist[-1]),
        "V_tf": float(V_hist[-1]),
    }

    if calcular_direto_H:
        dE_dH_m = float(_integral_trapezoidal(np.array(dP_dH_hist_m), dx=dt))

        dqin_adim_dH_m = (
            _row_dot_scalar(dA_scaled_sem_bc_dH_m, n_inlet, p_n)
            + _row_dot_scalar(A_scaled, n_inlet, sp_n)
        )
        dqin_dH_m = dqin_adim_dH_m * qref

        dV_dH_m = float(np.sum(sw_n * wref) * area_elemento_fisico)

        # Converte de derivada por metro para derivada por micrometro
        saida["dE_dH_direct"] = dE_dH_m * 1e-6
        saida["dqin_dH_direct"] = dqin_dH_m * 1e-6
        saida["dV_dH_direct"] = dV_dH_m * 1e-6

    return saida


def _derivada_forward(f0, fp, eps):
    return (fp - f0) / eps


def _derivada_centrada(fm, fp, eps):
    return (fp - fm) / (2.0 * eps)


def _sensibilidade_parametro(params_base, parametro, valores, eps, t_max, dt):
    """
    Calcula as derivadas por forward e centered difference para TC ou H.
    """
    E = []
    qin = []
    V = []

    dE_fwd = []
    dqin_fwd = []
    dV_fwd = []

    dE_ctr = []
    dqin_ctr = []
    dV_ctr = []

    dE_direct = []
    dqin_direct = []
    dV_direct = []

    for x in valores:
        print(f"   Rodando {parametro} = {x:.4g}")

        if parametro == "TC":
            base = SimulaCaso645(params_base, TC=x, H_um=params_base["H_k"] * 1e6,
                                 t_max=t_max, dt=dt)
            plus = SimulaCaso645(params_base, TC=x + eps, H_um=params_base["H_k"] * 1e6,
                                 t_max=t_max, dt=dt)
            minus = SimulaCaso645(params_base, TC=x - eps, H_um=params_base["H_k"] * 1e6,
                                  t_max=t_max, dt=dt)

            dE_direct.append(np.nan)
            dqin_direct.append(np.nan)
            dV_direct.append(np.nan)

        elif parametro == "H":
            base = SimulaCaso645(params_base, TC=params_base["TC"], H_um=x,
                                 t_max=t_max, dt=dt, calcular_direto_H=True)
            plus = SimulaCaso645(params_base, TC=params_base["TC"], H_um=x + eps,
                                 t_max=t_max, dt=dt)
            minus = SimulaCaso645(params_base, TC=params_base["TC"], H_um=x - eps,
                                  t_max=t_max, dt=dt)

            dE_direct.append(base["dE_dH_direct"])
            dqin_direct.append(base["dqin_dH_direct"])
            dV_direct.append(base["dV_dH_direct"])

        else:
            raise ValueError("parametro deve ser 'TC' ou 'H'.")

        E.append(base["E"])
        qin.append(base["qin_tf"])
        V.append(base["V_tf"])

        dE_fwd.append(_derivada_forward(base["E"], plus["E"], eps))
        dqin_fwd.append(_derivada_forward(base["qin_tf"], plus["qin_tf"], eps))
        dV_fwd.append(_derivada_forward(base["V_tf"], plus["V_tf"], eps))

        dE_ctr.append(_derivada_centrada(minus["E"], plus["E"], eps))
        dqin_ctr.append(_derivada_centrada(minus["qin_tf"], plus["qin_tf"], eps))
        dV_ctr.append(_derivada_centrada(minus["V_tf"], plus["V_tf"], eps))

    return {
        "parametro": parametro,
        "valores": np.array(valores, dtype=float),
        "E": np.array(E, dtype=float),
        "qin_tf": np.array(qin, dtype=float),
        "V_tf": np.array(V, dtype=float),
        "dE_forward": np.array(dE_fwd, dtype=float),
        "dqin_forward": np.array(dqin_fwd, dtype=float),
        "dV_forward": np.array(dV_fwd, dtype=float),
        "dE_centered": np.array(dE_ctr, dtype=float),
        "dqin_centered": np.array(dqin_ctr, dtype=float),
        "dV_centered": np.array(dV_ctr, dtype=float),
        "dE_direct": np.array(dE_direct, dtype=float),
        "dqin_direct": np.array(dqin_direct, dtype=float),
        "dV_direct": np.array(dV_direct, dtype=float),
    }


def RodaExercicio645_Sensibilidade(params_base, t_max=4.0, dt=0.05):
    """
    Resolve o Exercício 1 da Seção 6.4.5.

    Parâmetros:
        TC em [0, 250] °C
        H  em [500, 1500] micrometros

    Saídas:
        E, qin(tf), V(tf)
        derivadas forward
        derivadas centered
        derivadas diretas em relação a H
    """
    TC_grid = np.linspace(0.0, 250.0, 6)
    H_grid = np.linspace(500.0, 1500.0, 6)

    eps_TC = 1.0      # °C
    eps_H = 10.0      # micrometros

    print("\n[6.4.5] Sensibilidade em relação a TC")
    res_TC = _sensibilidade_parametro(
        params_base, "TC", TC_grid, eps_TC, t_max, dt
    )

    print("\n[6.4.5] Sensibilidade em relação a H")
    res_H = _sensibilidade_parametro(
        params_base, "H", H_grid, eps_H, t_max, dt
    )

    return {
        "TC": res_TC,
        "H": res_H,
        "eps_TC": eps_TC,
        "eps_H": eps_H,
        "t_max": t_max,
        "dt": dt,
    }


def PlotaExercicio645_Sensibilidade(resultados, save_dir=None):
    """
    Gera 5 figuras organizadas para o exercício 6.4.5:
      1) Saídas vs TC
      2) Sensibilidades vs TC
      3) Saídas vs H
      4) Sensibilidades vs H (forward e centered)
      5) Comparação centered vs método direto em H
    """

    # ==========================================================
    # FIGURA 1 - SAÍDAS EM FUNÇÃO DE TC
    # ==========================================================
    x_TC = resultados["TC"]["valores"]

    fig1, axs = plt.subplots(3, 1, figsize=(9, 11))

    axs[0].plot(x_TC, resultados["TC"]["E"], "o-")
    axs[0].set_title("Energia E em função de TC")
    axs[0].set_ylabel("E")
    axs[0].grid(True, linestyle=":", alpha=0.6)

    axs[1].plot(x_TC, resultados["TC"]["qin_tf"], "s-")
    axs[1].set_title(r"Vazão de entrada $q_{inlet}(t_f)$ em função de TC")
    axs[1].set_ylabel(r"$q_{inlet}(t_f)$")
    axs[1].grid(True, linestyle=":", alpha=0.6)

    axs[2].plot(x_TC, resultados["TC"]["V_tf"], "^-")
    axs[2].set_title(r"Volume final $V(t_f)$ em função de TC")
    axs[2].set_xlabel("TC [°C]")
    axs[2].set_ylabel(r"$V(t_f)$")
    axs[2].grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    if save_dir is not None:
        caminho = f"{save_dir}/ex645_fig1_saidas_vs_TC.png"
        plt.savefig(caminho, dpi=250, bbox_inches="tight")
        print(f"Figura salva: {caminho}")
    plt.show()

    # ==========================================================
    # FIGURA 2 - SENSIBILIDADES EM RELAÇÃO A TC
    # ==========================================================
    fig2, axs = plt.subplots(3, 1, figsize=(9, 11))

    axs[0].plot(x_TC, resultados["TC"]["dE_forward"], "o--", label="Forward")
    axs[0].plot(x_TC, resultados["TC"]["dE_centered"], "s-", label="Centered")
    axs[0].set_title(r"Sensibilidade $dE/dTC$")
    axs[0].set_ylabel(r"$dE/dTC$")
    axs[0].grid(True, linestyle=":", alpha=0.6)
    axs[0].legend()

    axs[1].plot(x_TC, resultados["TC"]["dqin_forward"], "o--", label="Forward")
    axs[1].plot(x_TC, resultados["TC"]["dqin_centered"], "s-", label="Centered")
    axs[1].set_title(r"Sensibilidade $dq_{inlet}/dTC$")
    axs[1].set_ylabel(r"$dq_{inlet}/dTC$")
    axs[1].grid(True, linestyle=":", alpha=0.6)
    axs[1].legend()

    axs[2].plot(x_TC, resultados["TC"]["dV_forward"], "o--", label="Forward")
    axs[2].plot(x_TC, resultados["TC"]["dV_centered"], "s-", label="Centered")
    axs[2].set_title(r"Sensibilidade $dV/dTC$")
    axs[2].set_xlabel("TC [°C]")
    axs[2].set_ylabel(r"$dV/dTC$")
    axs[2].grid(True, linestyle=":", alpha=0.6)
    axs[2].legend()

    plt.tight_layout()
    if save_dir is not None:
        caminho = f"{save_dir}/ex645_fig2_sensibilidade_vs_TC.png"
        plt.savefig(caminho, dpi=250, bbox_inches="tight")
        print(f"Figura salva: {caminho}")
    plt.show()

    # ==========================================================
    # FIGURA 3 - SAÍDAS EM FUNÇÃO DE H
    # ==========================================================
    x_H = resultados["H"]["valores"]

    fig3, axs = plt.subplots(3, 1, figsize=(9, 11))

    axs[0].plot(x_H, resultados["H"]["E"], "o-")
    axs[0].set_title("Energia E em função de H")
    axs[0].set_ylabel("E")
    axs[0].grid(True, linestyle=":", alpha=0.6)

    axs[1].plot(x_H, resultados["H"]["qin_tf"], "s-")
    axs[1].set_title(r"Vazão de entrada $q_{inlet}(t_f)$ em função de H")
    axs[1].set_ylabel(r"$q_{inlet}(t_f)$")
    axs[1].grid(True, linestyle=":", alpha=0.6)

    axs[2].plot(x_H, resultados["H"]["V_tf"], "^-")
    axs[2].set_title(r"Volume final $V(t_f)$ em função de H")
    axs[2].set_xlabel("H [µm]")
    axs[2].set_ylabel(r"$V(t_f)$")
    axs[2].grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    if save_dir is not None:
        caminho = f"{save_dir}/ex645_fig3_saidas_vs_H.png"
        plt.savefig(caminho, dpi=250, bbox_inches="tight")
        print(f"Figura salva: {caminho}")
    plt.show()

    # ==========================================================
    # FIGURA 4 - SENSIBILIDADES EM RELAÇÃO A H (FORWARD/CENTERED)
    # ==========================================================
    fig4, axs = plt.subplots(3, 1, figsize=(9, 11))

    axs[0].plot(x_H, resultados["H"]["dE_forward"], "o--", label="Forward")
    axs[0].plot(x_H, resultados["H"]["dE_centered"], "s-", label="Centered")
    axs[0].set_title(r"Sensibilidade $dE/dH$")
    axs[0].set_ylabel(r"$dE/dH$")
    axs[0].grid(True, linestyle=":", alpha=0.6)
    axs[0].legend()

    axs[1].plot(x_H, resultados["H"]["dqin_forward"], "o--", label="Forward")
    axs[1].plot(x_H, resultados["H"]["dqin_centered"], "s-", label="Centered")
    axs[1].set_title(r"Sensibilidade $dq_{inlet}/dH$")
    axs[1].set_ylabel(r"$dq_{inlet}/dH$")
    axs[1].grid(True, linestyle=":", alpha=0.6)
    axs[1].legend()

    axs[2].plot(x_H, resultados["H"]["dV_forward"], "o--", label="Forward")
    axs[2].plot(x_H, resultados["H"]["dV_centered"], "s-", label="Centered")
    axs[2].set_title(r"Sensibilidade $dV/dH$")
    axs[2].set_xlabel("H [µm]")
    axs[2].set_ylabel(r"$dV/dH$")
    axs[2].grid(True, linestyle=":", alpha=0.6)
    axs[2].legend()

    plt.tight_layout()
    if save_dir is not None:
        caminho = f"{save_dir}/ex645_fig4_sensibilidade_vs_H_fd.png"
        plt.savefig(caminho, dpi=250, bbox_inches="tight")
        print(f"Figura salva: {caminho}")
    plt.show()

    # ==========================================================
    # FIGURA 5 - COMPARAÇÃO CENTERED VS MÉTODO DIRETO EM H
    # ==========================================================
    fig5, axs = plt.subplots(3, 1, figsize=(9, 11))

    axs[0].plot(x_H, resultados["H"]["dE_centered"], "s-", label="Centered")
    axs[0].plot(x_H, resultados["H"]["dE_direct"], "^-", label="Direto")
    axs[0].set_title(r"Comparação $dE/dH$: Centered vs Método Direto")
    axs[0].set_ylabel(r"$dE/dH$")
    axs[0].grid(True, linestyle=":", alpha=0.6)
    axs[0].legend()

    axs[1].plot(x_H, resultados["H"]["dqin_centered"], "s-", label="Centered")
    axs[1].plot(x_H, resultados["H"]["dqin_direct"], "^-", label="Direto")
    axs[1].set_title(r"Comparação $dq_{inlet}/dH$: Centered vs Método Direto")
    axs[1].set_ylabel(r"$dq_{inlet}/dH$")
    axs[1].grid(True, linestyle=":", alpha=0.6)
    axs[1].legend()

    axs[2].plot(x_H, resultados["H"]["dV_centered"], "s-", label="Centered")
    axs[2].plot(x_H, resultados["H"]["dV_direct"], "^-", label="Direto")
    axs[2].set_title(r"Comparação $dV/dH$: Centered vs Método Direto")
    axs[2].set_xlabel("H [µm]")
    axs[2].set_ylabel(r"$dV/dH$")
    axs[2].grid(True, linestyle=":", alpha=0.6)
    axs[2].legend()

    plt.tight_layout()
    if save_dir is not None:
        caminho = f"{save_dir}/ex645_fig5_comparacao_metodo_direto.png"
        plt.savefig(caminho, dpi=250, bbox_inches="tight")
        print(f"Figura salva: {caminho}")
    plt.show()

# ========================================================================
# SEÇÃO 6.4.5 - EXERCÍCIO 2
#  Localização de Raízes via Newton-Raphson
# ========================================================================

import numpy as np
import matplotlib.pyplot as plt

# Nota: Certifique-se de que a função SimulaCaso645 está importada ou definida neste arquivo
# ou que o seu código original já esteja lidando com ela!

def funcao_energia_H(GD_params, H_um, E_alvo=7.5):
    """
    Calcula f(H) = E(H) - E_alvo.
    """
    resultado = SimulaCaso645(GD_params, H_um=H_um)
    return resultado["E"] - E_alvo

def derivada_forward_energia(GD_params, H_um, E_alvo=7.5, eps=1.0):
    """
    Aproxima f'(H) por diferença finita progressiva.
    """
    f_H = funcao_energia_H(GD_params, H_um, E_alvo)
    f_H_plus = funcao_energia_H(GD_params, H_um + eps, E_alvo)

    return (f_H_plus - f_H) / eps

def newton_raphson_energia(
    GD_params,
    H_chute,
    E_alvo=7.5,
    tol=1e-4,
    max_iter=10,
    eps=1.0,
    H_min=500.0,
    H_max=1500.0
):
    """
    Resolve f(H) = E(H) - E_alvo = 0 pelo método de Newton-Raphson.
    Retorna o H encontrado e o histórico das iterações.
    """

    H = float(H_chute)
    historico = []

    for i in range(max_iter):
        f_H = funcao_energia_H(GD_params, H, E_alvo)
        df_H = derivada_forward_energia(GD_params, H, E_alvo, eps)

        historico.append({
            "iter": i,
            "H": H,
            "f_H": f_H,
            "df_H": df_H
        })

        if abs(f_H) < tol:
            return H, historico

        if abs(df_H) < 1e-12:
            raise RuntimeError("Derivada muito próxima de zero. Newton-Raphson falhou.")

        H_novo = H - f_H / df_H

        # Evita sair do intervalo físico do problema
        H_novo = max(H_min, min(H_max, H_novo))

        H = H_novo

    return H, historico

def plotar_raiz_energia(GD_params, H_otimo, E_alvo=7.5):
    """
    Plota f(H) = E(H) - E_alvo no intervalo físico considerado.
    """

    Hs = np.linspace(500, 1500, 50)
    fs = []

    for H in Hs:
        f_H = funcao_energia_H(GD_params, H, E_alvo)
        fs.append(f_H)

    plt.figure(figsize=(8, 5))
    plt.plot(Hs, fs, label="f(H) = E(H) - 7.5")
    plt.axhline(0, linestyle="--", color="black")
    plt.axvline(H_otimo, linestyle="--", label=f"H = {H_otimo:.2f} µm")
    plt.scatter([H_otimo], [0], zorder=5)

    plt.xlabel("H [µm]")
    plt.ylabel("f(H)")
    plt.title("Localização da raiz por Newton-Raphson")
    plt.grid(True)
    plt.legend()
    plt.show()
    plt.legend()
    plt.savefig("grafico_raiz_newton.png")
# =======================================================================================
# CONCLUSÃO
# =======================================================================================

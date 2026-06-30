"""
functionsGD.py
----------------------------------------------------------------------
Funcoes de montagem da BASE FISICA do GEMEO DIGITAL (Capitulo 6).

Este modulo NAO resolve os exercicios das Secoes 6.3 (Predicoes sob
Incerteza - Amostragem Monte Carlo) e 6.4 (Aprendizado de Modelos -
Interpolacao, Minimos Quadrados, Zeros de Funcao e Diferenciacao
Numerica). Ele concentra apenas as rotinas necessarias para montar, de
ponta a ponta, o modelo acoplado completo (Placa Termica -> Rede
Hidraulica -> Membrana Elastica), nas condicoes nominais especificadas
na Secao 6.2 do PDF da disciplina.

Pipeline fisico (conforme Secao 6.1 do PDF):

    1) Resolve-se o balanco termico da placa solida, considerando que a
       proximidade dos microcanais modifica localmente a condutividade
       termica efetiva do material (Equacao 4.3, ja implementada em
       ObterCondutividadeFaces_ViaNos no Capitulo 4);

    2) Estima-se a temperatura media de cada canal (integral de linha
       do campo de temperatura ao longo de cada aresta da rede) e, a
       partir dela, a viscosidade local do fluido e a condutancia
       hidraulica atualizada de cada canal;

    3) Monta-se a matriz de rigidez hidraulica da rede com as
       condutancias termicamente atualizadas e acopla-se, de forma
       monolitica (conforme estruturado no Capitulo 5), a membrana
       elastica, resultando na matriz global do Gemeo Digital.

NOTA SOBRE OS IMPORTS: seguindo o mesmo padrao ja usado em
functionsHM.py (imports "achatados", sem prefixo de pacote), este
arquivo importa diretamente "functions", "functionsM", "functionsHT" e
"functionsHM". Para isso funcionar, o script que importar functionsGD
(o mainGD.py) precisa ter adicionado as pastas RedeHidraulica,
MembranaElastica, AcoplamentoHidraulicoTermico e
AcoplamentoHidraulicoMecanico diretamente ao sys.path antes do import
- o que ja e feito em mainGD.py.
"""

import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Imports cruzados entre os modulos dos capitulos anteriores (imports
# "achatados", no mesmo padrao usado internamente por functionsHM.py)
# ----------------------------------------------------------------------
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


# ========================================================================
# 7. SECAO 6.3.2 - EXERCICIO 1
#    Analise Estacionaria de Falhas Hidraulicas (Amostragem Monte Carlo)
# ========================================================================
#
# Cenario fisico (Secao 6.3 / 6.3.2 do PDF):
#   - Subsistema ISOLADO composto apenas pela Rede Hidraulica e pela
#     Placa Termica, em regime permanente (sem a membrana elastica).
#   - Cada canal possui probabilidade independente p_O de obstrucao;
#     quando obstruido, sua condutancia C_k e reduzida por um fator de
#     severidade f_obs:  C_k <- C_k / f_obs.
#   - Pressao de descarga nula no reservatorio: p_outlet = 0.
#   - Pergunta de projeto: Prob(q_inlet < q_critico),
#     com q_critico = 1.25e-5 (vazao volumetrica no no 0).
#
# ------------------------------------------------------------------------
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
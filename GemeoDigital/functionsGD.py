import numpy as np

from RedeHidraulica import functions as fH
from PlacaTermica import functionsT as fT
from MembranaElastica import functionsM as fM

def RandomFail(condutancias_nominais, p0, fobs):
    """
    Aplica falhas locais randômicas nas condutâncias da rede microfluídica.
    
    Parâmetros:
    - condutancias_nominais: array com as condutâncias (C) de todos os canais da rede
    - p0: probabilidade (entre 0 e 1) de um canal sofrer obstrução
    - fobs: fator de obstrução, o quão severa é a falha (ex: 0.5 reduz a condutância pela metade)
    """
    C_atualizado = condutancias_nominais.copy()
    n_canais = len(C_atualizado)
    
    for i in range(n_canais):
        # Sorteia um número entre 0 e 1 para cada canal
        if np.random.rand() < p0:
            C_atualizado[i] = C_atualizado[i] * fobs
            
    return C_atualizado


def ResolverGemeoDigital(params_placa, params_rede, params_membrana, falha_kwargs=None):
    """
    Resolve a cadeia multifísica do Gêmeo Digital:
    Passo 1 (Placa) -> Passo 2 (Rede Hidráulica) -> Passo 3 (Membrana)
    
    Retorna os principais Outputs/QoI (Quantities of Interest).
    """
    
    # ---------------------------------------------------------
    # PASSO 1: PLACA TÉRMICA
    # ---------------------------------------------------------
    # Resolve a placa para encontrar o campo de temperatura
    # Adapte os nomes das funções (ex: resolver_placa_2d) para os que você criou no Cap 4
    T_campo, malha_X, malha_Y = fT.resolver_placa(**params_placa)
    
    # Mapeamento do campo de temperatura para o centro de cada microcanal
    # (Ou temperatura média da placa afetando os canais, conforme a sua implementação do Cap 4/6)
    T_canais = fT.extrair_temperatura_nos_canais(T_campo, malha_X, malha_Y, params_rede['coords_canais'])
    
    
    # ---------------------------------------------------------
    # PASSO 2: REDE HIDRÁULICA E ACOPLAMENTO TÉRMICO
    # ---------------------------------------------------------
    # Recalcula a viscosidade para cada canal usando a temperatura extraída
    viscosidades = fH.calcular_viscosidade_fluido(T_canais)
    
    # Determina as condutâncias base da rede com a viscosidade atualizada
    C_canais = fH.calcular_condutancias(viscosidades, params_rede['largura_canal'])
    
    # Condição para os testes de Monte Carlo: aplicar a falha antes da resolução
    if falha_kwargs is not None:
        C_canais = RandomFail(C_canais, falha_kwargs['p0'], falha_kwargs['fobs'])
        
    # Resolve as pressões dos nós e vazões
    p_nos, q_canais = fH.resolver_rede_pressao(
        C_canais, 
        params_rede['matriz_incidencia'], 
        params_rede['p_inlet'], 
        params_rede['no_inlet'], 
        params_rede['no_outlet']
    )
    
    # Extrai QoI (Vazão total no inlet e Pressão no outlet)
    # Supondo que a vazão do inlet seja a soma das vazões dos canais conectados ao nó 0
    q_inlet_total = np.sum(q_canais[params_rede['indices_canais_inlet']])
    p_outlet = p_nos[params_rede['no_outlet']]
    
    
    # ---------------------------------------------------------
    # PASSO 3: MEMBRANA ELÁSTICA E ACOPLAMENTO MECÂNICO
    # ---------------------------------------------------------
    # A pressão de outlet da rede entra como carga na membrana
    params_membrana['pressao_aplicada'] = p_outlet
    
    # Resolve a deflexão transiente no tempo
    tempo, W_historico = fM.resolver_membrana_transiente(**params_membrana)
    
    # Extrai o QoI final: Deflexão central (w no nó 0) ao longo do tempo ou seu valor máximo/estacionário
    w_centro_historico = W_historico[:, 0]
    
    return {
        'T_campo': T_campo,
        'q_inlet_total': q_inlet_total,
        'p_outlet': p_outlet,
        'tempo': tempo,
        'w_centro': w_centro_historico,
        'w_max_absoluto': np.max(np.abs(w_centro_historico))
    }
import numpy as np
from scipy import sparse
import matplotlib.pyplot as plt

def ij2n(i, j, N1):
    """Converte coordenadas 2D (i, j) para índice 1D."""
    return i + j * N1

def BuildMatrizes_Eigen(N1, N2, sigma, rho, e, h):
    """
    Monta as matrizes de Rigidez (K) e Massa (M) para o problema de autovalores
    generalizado da membrana elástica retangular.
    """
    nunk = N1 * N2
    
    # 1. Matriz de Rigidez (K)
    d0 = 4.0 * np.ones(nunk)
    d1 = -np.ones(nunk - 1)
    dN = -np.ones(nunk - N1)
    
    # Quebra de conexão entre linhas para evitar wrap-around
    for i in range(1, N2):
        d1[i * N1 - 1] = 0
        
    # Inicializando K em formato LIL (List of Lists) para facilitar a modificação
    K = (sigma / h**2) * sparse.diags(
        [dN, d1, d0, d1, dN], 
        [-N1, -1, 0, 1, N1], 
        format='lil'
    )
    
    # Penalização para forçar deslocamento nulo nas bordas (Dirichlet)
    big_number = 10000.0
    
    # Identificação dos nós de contorno
    nos_contorno = []
    for i in range(N1):
        for j in range(N2):
            # Se for borda (esquerda, direita, base, topo)
            if i == 0 or i == N1 - 1 or j == 0 or j == N2 - 1:
                nos_contorno.append(ij2n(i, j, N1))
                
    # Aplicação da penalidade (cravando big_number na diagonal e isolando o nó)
    for Ic in nos_contorno:
        K[Ic, :] = 0
        K[:, Ic] = 0
        K[Ic, Ic] = big_number
        
    # Converte K para CSR para cálculos eficientes com álgebra linear esparsa
    K = K.tocsr()
    
    # 2. Matriz de Massa (M)
    # Como a massa é uniformemente distribuída, M é múltipla da matriz Identidade
    M = rho * e * sparse.identity(nunk, format='csr')
    
    return K, M

# Exercício 1
def BuildMatrizes_Eigen_Circular(N1, N2, sigma_ad, rho_ad, e_ad, h):
    """
    Constrói K e M para a membrana circular adimensionalizada.
    O domínio físico simulado é um quadrado 2x2.
    """
    nunk = N1 * N2
    
    # 1. Matriz de Rigidez K (Laplaciano)
    d0 = 4.0 * np.ones(nunk)
    d1 = -np.ones(nunk - 1)
    dN = -np.ones(nunk - N1)
    
    # Corte estrutural do wrap-around
    for i in range(1, N2):
        d1[i * N1 - 1] = 0
        
    K = (sigma_ad / h**2) * sparse.diags(
        [dN, d1, d0, d1, dN], 
        [-N1, -1, 0, 1, N1], 
        format='lil'
    )
    
    # 2. Matriz de Massa M
    M = rho_ad * e_ad * sparse.identity(nunk, format='csr')
    
    # 3. Adimensionalização e Geração da Máscara Circular
    Lx = 2.0
    Ly = 2.0
    R_ad = 1.0
    xc = 1.0
    yc = 1.0
    
    x = np.linspace(0, Lx, N1)
    y = np.linspace(0, Ly, N2)
    
    big_number = 1e4
    nos_restritos = []
    
    # Identificação dos nós exteriores à membrana
    for i in range(N1):
        for j in range(N2):
            dist_quadrada = (x[i] - xc)**2 + (y[j] - yc)**2
            
            # Se o nó estiver fora do círculo geométrico ou exatamente na borda
            if dist_quadrada >= R_ad**2:
                nos_restritos.append(ij2n(i, j, N1))
                
    # Injeção da condição de Dirichlet com penalização pesada
    for Ic in nos_restritos:
        K[Ic, :] = 0
        K[:, Ic] = 0
        K[Ic, Ic] = big_number
        
    return K.tocsr(), M

def PlotaModo(N1, N2, Lx, Ly, phi, modo_idx, omega):
    """
    Plota a superfície 3D de um modo de vibração específico da membrana.
    """
    x = np.linspace(0, Lx, N1)
    y = np.linspace(0, Ly, N2)
    X, Y = np.meshgrid(x, y)
    
    # Remodela o autovetor 1D de volta para o grid 2D
    Z = phi.reshape((N2, N1))
    
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')
    
    ax.set_title(f'Modo Fundamental {modo_idx} \n Frequência $\\omega$ = {omega:.2f} rad/s')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Deslocamento (w)')
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5, label='Amplitude Adimensional')
    
    plt.tight_layout()
    plt.show()

# Exercício 4

def DecomposicaoModal(Phi, M, Z):
    
    nmodes = Phi.shape[1]
    alpha = np.zeros(nmodes)

    for i in range(nmodes):

        phi_i = Phi[:, i]

        numerador = phi_i.T @ Z
        denominador = phi_i.T @ (M @ phi_i)

        alpha[i] = numerador / denominador

    return alpha
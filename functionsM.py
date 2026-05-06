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
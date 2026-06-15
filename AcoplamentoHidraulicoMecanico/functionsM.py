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
    
    # 1. Matrizes Base em CSR
    d0 = 4.0 * np.ones(nunk)
    d1 = -np.ones(nunk - 1)
    dN = -np.ones(nunk - N1)
    
    K = (sigma / h**2) * sparse.diags(
        [dN, d1, d0, d1, dN], 
        [-N1, -1, 0, 1, N1], 
        format='csr'
    )
    M = rho * e * sparse.identity(nunk, format='csr')
    
    # 2. Configuração da Penalidade
    big_number = 10000.0
    Iden = big_number * sparse.identity(nunk, format='csr')
    
    # 3. Aplicação de Dirichlet por Slicing (Lados Horizontais e Verticais)
    # Lados verticais (Esquerda e Direita)
    for j in range(0, N2):
        # Esquerda (i=0)
        Ic_L = ij2n(0, j, N1)
        K[Ic_L, :], K[:, Ic_L] = Iden[Ic_L, :], Iden[:, Ic_L]
        
        # Direita (i=N1-1)
        Ic_R = ij2n(N1 - 1, j, N1)
        K[Ic_R, :], K[:, Ic_R] = Iden[Ic_R, :], Iden[:, Ic_R]
        
    # Lados horizontais (Base e Topo)
    for i in range(0, N1):
        # Base (j=0)
        Ic_B = ij2n(i, 0, N1)
        K[Ic_B, :], K[:, Ic_B] = Iden[Ic_B, :], Iden[:, Ic_B]
        
        # Topo (j=N2-1)
        Ic_T = ij2n(i, N2 - 1, N1)
        K[Ic_T, :], K[:, Ic_T] = Iden[Ic_T, :], Iden[:, Ic_T]
        
    return K, M

# Exercício 1
def BuildMatrizes_Eigen_Circular(N1, N2, sigma_ad, rho_ad, e_ad, h):
    """
    Constrói K e M para a membrana circular adimensionalizada.
    O domínio físico simulado é um quadrado 2x2.
    """
 
    nunk = N1 * N2
       
    # 1. Matrizes Base em CSR
    d0 = 4.0 * np.ones(nunk)
    d1 = -np.ones(nunk - 1)
    dN = -np.ones(nunk - N1)
    
    K = (sigma_ad / h**2) * sparse.diags(
        [dN, d1, d0, d1, dN], 
        [-N1, -1, 0, 1, N1], 
        format='csr'
    )
    M = rho_ad * e_ad * sparse.identity(nunk, format='csr')
    
    # 2. Penalidade via Slicing
    big_number = 10000.0
    Iden = big_number * sparse.identity(nunk, format='csr')
    
    Lx, Ly = 2.0, 2.0
    R_ad, xc, yc = 1.0, 1.0, 1.0
    x = np.linspace(0, Lx, N1)
    y = np.linspace(0, Ly, N2)
    
    for j in range(N2):
        for i in range(N1):
            dist_quadrada = (x[i] - xc)**2 + (y[j] - yc)**2
            if dist_quadrada >= R_ad**2:
                Ic = ij2n(i, j, N1)
                # Atribuição direta de linha e coluna
                K[Ic, :], K[:, Ic] = Iden[Ic, :], Iden[:, Ic]
                
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
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union
from scipy import sparse
from scipy.sparse.linalg import spsolve
import time

def ij2n (i, j, Nx):
    return i + j*Nx

# FUNÇÃO ASSEMBLY --------------------------------------------------------------------------
import numpy as np

def Assembly(Nx, Ny, k):
    
    nunk = Nx * Ny
    A = np.zeros((nunk, nunk))
    
    # Percorre todos os pontos
    for i in range(Nx):
        for j in range(Ny):
            
            Ic = ij2n(i, j, Nx)
            
            # 🔹 Condição de contorno
            if i == 0 or i == Nx-1 or j == 0 or j == Ny-1:
                A[Ic, Ic] = 1
            
            else:
                # vizinhos
                Ie = ij2n(i+1, j, Nx)
                Iw = ij2n(i-1, j, Nx)
                In = ij2n(i, j+1, Nx)
                Is = ij2n(i, j-1, Nx)
                
                A[Ic, Ic] = 4*k
                A[Ic, Ie] = -k
                A[Ic, Iw] = -k
                A[Ic, In] = -k
                A[Ic, Is] = -k

    return A

# FUNÇÃO SolveSystem --------------------------------------------------------------------------
import numpy as np

def SolveSystem(Nx, Ny, h, k, TL, TR, TB, TT, fonte):
    
    # Tempo do Assembly
    t0 = time.time()
    A = Assembly(Nx, Ny, k)
    t_assembly = time.time() - t0

    # Tempo de montagem do sistema
    t0 = time.time()

    Atilde = A.copy()
    nunk = Nx * Ny
    b = np.zeros(nunk)
    
    if fonte is not None:
        b += fonte * h**2
    
    # 3. Matriz identidade (para impor contorno)
    Iden = np.identity(nunk)
    
    # 4. Aplicar condições de contorno
    for i in range(Nx):
        for j in range(Ny):
            
            Ic = ij2n(i, j, Nx)
            
            # esquerda
            if i == 0:
                Atilde[Ic, :] = Iden[Ic, :]
                b[Ic] = TL
            
            # direita
            elif i == Nx-1:
                Atilde[Ic, :] = Iden[Ic, :]
                b[Ic] = TR
            
            # base
            elif j == 0:
                Atilde[Ic, :] = Iden[Ic, :]
                b[Ic] = TB[i]
            
            # topo
            elif j == Ny-1:
                Atilde[Ic, :] = Iden[Ic, :]
                b[Ic] = TT[i]

    t_montagem = time.time() - t0
    
    # tempo reolução do sistema linear
    t0 = time.time()
    T = np.linalg.solve(Atilde, b)
    t_sistema = time.time() - t0
    
    T_grid = T.reshape((Ny, Nx))

    return T_grid, t_assembly, t_montagem, t_sistema

# FUNÇÃO SolveSystemSparse --------------------------------------------------------------------------
def SolveSystemSparse(Nx, Ny, h, k, TL, TR, TB, TT, fonte):


    nunk = Nx * Ny
    
    #tempo para a montagem da matriz
    t0 = time.time()
    # 1. Diagonais
    d0 = np.ones(nunk) * 4.0 * k          # diagonal principal
    d1 = -np.ones(nunk - 1) * k           # vizinhos esquerda/direita
    dN = -np.ones(nunk - Nx) * k          # vizinhos cima/baixo
    
    for i in range(1, Ny):
        d1[i * Nx - 1] = 0  # quebra conexão entre linhas
    
    # 2. Montar matriz esparsa
    A = sparse.diags(
        [dN, d1, d0, d1, dN],
        [-Nx, -1, 0, 1, Nx],
        format='lil'  # fácil de modificar
    )

    t_assembly = time.time() - t0
    
    # tempo para a montagem do sistema:
    t0 = time.time()

    # 3. Vetor do lado direito
    b = np.zeros(nunk)
    
    if fonte is not None:
        b += fonte * h**2
    
    # 4. Aplicar condições de contorno
    for i in range(Nx):
        for j in range(Ny):
            
            Ic = ij2n(i, j, Nx)
            
            if i == 0:  # esquerda
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TL
                
            elif i == Nx-1:  # direita
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TR
                
            elif j == 0:  # base
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TB[i]
                
            elif j == Ny-1:  # topo
                A[Ic, :] = 0
                A[Ic, Ic] = 1
                b[Ic] = TT[i]
    
    # 5. Converter para CSR (eficiente)
    A = A.tocsr()

    t_montagem = time.time() - t0

    # tempo reolução do sistema linear
    t0 = time.time()
    
    # 6. Resolver sistema
    T = spsolve(A, b)

    t_sistema = time.time() - t0
    
    # 7. Converter para matriz 2D
    T_grid = T.reshape((Ny, Nx))
    
    return T_grid, t_assembly, t_montagem, t_sistema


# FUNÇÃO PlotaPlaca --------------------------------------------------------------------------

def PlotaPlaca(Nx, Ny, Lx, Ly, T, flag_type='contour', filename=None):
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
    Z = np.copy(T).reshape(Ny, Nx)
    if(flag_type == 'contour'):
      fig, ax = plt.subplots(figsize=(6,6))
      ax.set_aspect('equal')
      ax.set(xlabel='x', ylabel='y', title='Contours of temperature')
      im = ax.contourf(X, Y, Z, 20, cmap='jet')
      im2 = ax.contour(X, Y, Z, 20, linewidths=0.25, colors='k')
      fig.colorbar(im, ax=ax, orientation='horizontal')
    elif(flag_type == 'surface'):
      fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
      ax.set_aspect('equal')
      surf = ax.plot_surface(X, Y, Z, cmap='jet')
      fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5) 
    
    plt.xticks([0, Lx/2, Lx])
    plt.yticks([0, Ly/2, Ly])

    if(filename is not None):
      plt.savefig(filename)

    plt.show()

    return

#FUNÇÃO SolveSystemSparse_Circle ---------------------------------------------------------------------------------------
def SolveSystemSparse_Circle(Nx, Ny, h, k, TL, TR, TB, TT, fonte, Lx, Ly, R, xc, yc, TC):
 
    nunk = Nx * Ny
 
    # Coordenadas dos pontos da malha
    x_coords = np.linspace(0, Lx, Nx)
    y_coords = np.linspace(0, Ly, Ny)
 
    circle_mask = np.zeros((Ny, Nx), dtype=bool)
    for j in range(Ny):
        for i in range(Nx):
            dist = np.sqrt((x_coords[i] - xc)**2 + (y_coords[j] - yc)**2)
            if dist <= R:
                circle_mask[j, i] = True
 
    # Tempo de Assembly 
    t0 = time.time()
 
    d0 = np.ones(nunk)     *  4.0 * k
    d1 = np.ones(nunk - 1) * -k
    dN = np.ones(nunk - Nx) * -k
 
    for i in range(1, Ny):
        d1[i * Nx - 1] = 0  
 
    A = sparse.diags(
        [dN, d1, d0, d1, dN],
        [-Nx, -1, 0, 1, Nx],
        format='lil'
    )
 
    t_assembly = time.time() - t0
 
    # Tempo de montagem do sistema 
    t0 = time.time()
 
    b = np.zeros(nunk)


    for j in range(Ny):
        for i in range(Nx):
            Ic = ij2n(i, j, Nx)
            is_border = (i == 0 or i == Nx - 1 or j == 0 or j == Ny - 1)
            if fonte is not None and not is_border and not circle_mask[j, i]:
                b[Ic] = fonte * h ** 2

    for i in range(Nx):
        for j in range(Ny):
 
            Ic = ij2n(i, j, Nx)
 
            if i == 0:           # borda esquerda
                A[Ic, :] = 0;  A[Ic, Ic] = 1;  b[Ic] = TL
 
            elif i == Nx - 1:    # borda direita
                A[Ic, :] = 0;  A[Ic, Ic] = 1;  b[Ic] = TR
 
            elif j == 0:         # borda inferior
                A[Ic, :] = 0;  A[Ic, Ic] = 1;  b[Ic] = TB[i]
 
            elif j == Ny - 1:    # borda superior
                A[Ic, :] = 0;  A[Ic, Ic] = 1;  b[Ic] = TT[i]
 
            elif circle_mask[j, i]:   # região circular interna
                A[Ic, :] = 0;  A[Ic, Ic] = 1;  b[Ic] = TC
 
    A = A.tocsr()
    t_montagem = time.time() - t0
 
    # Tempo de resolução
    t0 = time.time()
    T  = spsolve(A, b)
    t_sistema = time.time() - t0
 
    T_grid = T.reshape((Ny, Nx))
 
    return T_grid, t_assembly, t_montagem, t_sistema, circle_mask

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import cm
import pandas as pd
from shapely.geometry import LineString, Point
from shapely.ops import unary_union
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import solve_triangular
from scipy.sparse.linalg import splu
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

    # 4. Aplicar condições de contorno (Abordagem Simétrica)   

    nos_contorno = []
    valores_contorno = []

    for i in range(Nx):
        for j in range(Ny):
            Ic = ij2n(i, j, Nx)
            if i == 0:
                nos_contorno.append(Ic); valores_contorno.append(TL)
            elif i == Nx - 1:
                nos_contorno.append(Ic); valores_contorno.append(TR)
            elif j == 0:
                nos_contorno.append(Ic); valores_contorno.append(TB[i])
            elif j == Ny - 1:
                nos_contorno.append(Ic); valores_contorno.append(TT[i])

    nos_contorno = np.array(nos_contorno)
    valores_contorno = np.array(valores_contorno)

    A_csc = A.tocsc()

    for idx, Ic in enumerate(nos_contorno):
        b -= A_csc[:, Ic].toarray().flatten() * valores_contorno[idx]

    A = A_csc.tolil()
    
    for Ic in nos_contorno:
        A[Ic, :] = 0  # Zera a linha
        A[:, Ic] = 0  # Zera a coluna (Garante a Simetria!)

    for idx, Ic in enumerate(nos_contorno):
        A[Ic, Ic] = 1
        b[Ic] = valores_contorno[idx]

    # 5. Converter para CSR (eficiente)
    A = A.tocsr() # Matriz simétrica final pronta para resolução

    t_montagem = time.time() - t0

    # tempo reolução do sistema linear
    t0 = time.time()
    
    # 6. Resolver sistema
    T = spsolve(A, b) # a função spsolve identifica que A é simétrica e usa um método otimizado para resolver o sistema

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
      ax.set(xlabel='x', ylabel='y', title=f'Contours of temperature ({Nx}x{Ny})')
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

    nos_fixos = []
    valores_fixos = []

    for j in range(Ny):
        for i in range(Nx):
            Ic = ij2n(i, j, Nx)
            
            # Avalia primeiro se é contorno ou círculo
            if i == 0:
                nos_fixos.append(Ic); valores_fixos.append(TL)
            elif i == Nx - 1:
                nos_fixos.append(Ic); valores_fixos.append(TR)
            elif j == 0:
                nos_fixos.append(Ic); valores_fixos.append(TB[i])
            elif j == Ny - 1:
                nos_fixos.append(Ic); valores_fixos.append(TT[i])
            elif circle_mask[j, i]:
                nos_fixos.append(Ic); valores_fixos.append(TC)
            
            # Se não é contorno nem círculo, é nó livre (aplica a fonte)
            else:
                if fonte is not None:
                    b[Ic] = fonte * h ** 2

    nos_fixos = np.array(nos_fixos)
    valores_fixos = np.array(valores_fixos)

    A_csc = A.tocsc()
    for idx, Ic in enumerate(nos_fixos):
        b -= A_csc[:, Ic].toarray().flatten() * valores_fixos[idx]

    A = A_csc.tolil()
    
    for Ic in nos_fixos:
        A[Ic, :] = 0
        A[:, Ic] = 0
        
    for idx, Ic in enumerate(nos_fixos):
        A[Ic, Ic] = 1 
        b[Ic] = valores_fixos[idx]

    A = A.tocsr()
    t_montagem = time.time() - t0
 
    # Tempo de resolução
    t0 = time.time()
    T  = spsolve(A, b)
    t_sistema = time.time() - t0
 
    T_grid = T.reshape((Ny, Nx))
 
    return T_grid, t_assembly, t_montagem, t_sistema, circle_mask

# FUNÇÃO Prepara_Sistema_Otimizado ---------------------------------------------------------
# Essa monta a matriz e fatora (LU) uma única vez, para ser usada em métodos iterativos onde a matriz não muda.

def Prepara_Sistema_Otimizado(Nx, Ny, h, k, TL, TR, TB, TT, fonte, Lx, Ly, R, xc, yc):

    nunk = Nx * Ny
    x_coords = np.linspace(0, Lx, Nx)
    y_coords = np.linspace(0, Ly, Ny)

    circle_mask = np.zeros((Ny, Nx), dtype=bool)
    for j in range(Ny):
        for i in range(Nx):
            dist = np.sqrt((x_coords[i] - xc)**2 + (y_coords[j] - yc)**2)
            if dist <= R:
                circle_mask[j, i] = True
    
    # Montagem da matriz A (Simétrica e Esparsa)
    d0, d1, dN = np.ones(nunk)*4*k, -np.ones(nunk-1)*k, -np.ones(nunk-Nx)*k
    for i in range(1, Ny): d1[i*Nx-1] = 0
    A = sparse.diags([dN, d1, d0, d1, dN], [-Nx, -1, 0, 1, Nx], format='lil')
    
    # Vetor b_base (Fonte de calor e contornos fixos, SEM o círculo ainda)
    b_base = np.zeros(nunk)
    nos_fixos = []
    
    for j in range(Ny):
        for i in range(Nx):
            Ic = ij2n(i, j, Nx)
            if i == 0: nos_fixos.append(Ic); b_base[Ic] = TL
            elif i == Nx-1: nos_fixos.append(Ic); b_base[Ic] = TR
            elif j == 0: nos_fixos.append(Ic); b_base[Ic] = TB[i]
            elif j == Ny-1: nos_fixos.append(Ic); b_base[Ic] = TT[i]
            elif circle_mask[j, i]: nos_fixos.append(Ic)
            else: b_base[Ic] = fonte * h**2

    # Aplicar zeros na matriz para os nós fixos (Dirichlet)
    for Ic in nos_fixos:
        A[Ic, :] = 0; A[Ic, Ic] = 1
        
    A_fatorada = splu(A.tocsc()) 
    
    return A_fatorada, b_base, circle_mask.flatten()

# FUNÇÃO Resolve_Rapido --------------------------------------------------------------------------
# Resolve o sistemae em milissegundos usando a matriz já fatorada. Ideal para métodos iterativos onde só o lado direito muda (ex: TC do círculo).

def Resolve_Rapido(A_fatorada, b_base, circle_mask_1d, TC_atual):

    b = b_base.copy()
    b[circle_mask_1d] = TC_atual # Atualiza apenas os pontos do círculo
    
    T = A_fatorada.solve(b) # Resolve por substituição
    return T

# FUNÇÃO Jacobi --------------------------------------------------------------------------

def Jacobi(Nx, Ny, h, k, TL, TR, TB, TT, fonte, TOL, MAXIT, animation, frame_skip):
    nunk = Nx * Ny
    
    # 1. Montagem da Matriz Esparsa (Idêntico ao SolveSystemSparse)
    d0 = np.ones(nunk) * 4.0 * k
    d1 = -np.ones(nunk - 1) * k
    dN = -np.ones(nunk - Nx) * k
    for i in range(1, Ny): d1[i * Nx - 1] = 0
    A = sparse.diags([dN, d1, d0, d1, dN], [-Nx, -1, 0, 1, Nx], format='lil')
    
    b = np.zeros(nunk)
    
    # 2. Único laço para Contornos e Fonte
    for i in range(Nx):
        for j in range(Ny):
            Ic = ij2n(i, j, Nx)
            if i == 0:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TL
            elif i == Nx-1:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TR
            elif j == 0:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TB[i]
            elif j == Ny-1:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TT[i]
            else:
                if fonte is not None:
                    b[Ic] = fonte * h**2
                    
    A = A.tocsr()
    
    # 3. Preparação do Jacobi
    diag_A = A.diagonal() 
    x = np.zeros(nunk)
    F = A.dot(x) - b   
    
    k_iter = 0
    frames = [] if animation else None
    t0 = time.time()
    
    # 4. Loop Iterativo
    while np.linalg.norm(F, np.inf) > TOL and k_iter < MAXIT:
        d = -F / diag_A   # Divisão Vetorial Direta
        x = x + d
        F = A.dot(x) - b
        k_iter += 1
        
        if animation and (k_iter % frame_skip == 0):
            frames.append(x.reshape((Ny, Nx)).copy())

    tempo = time.time() - t0
    T_grid = x.reshape((Ny, Nx))
    return T_grid, k_iter, tempo, frames

# FUNÇÃO Gauss-Seidel --------------------------------------------------------------------------

from scipy.sparse.linalg import spsolve_triangular

def GaussSeidel(Nx, Ny, h, k, TL, TR, TB, TT, fonte, TOL, MAXIT, animation, frame_skip):
    nunk = Nx * Ny
    
    # 1. Montagem da Matriz Esparsa
    d0 = np.ones(nunk) * 4.0 * k
    d1 = -np.ones(nunk - 1) * k
    dN = -np.ones(nunk - Nx) * k
    for i in range(1, Ny): d1[i * Nx - 1] = 0
    A = sparse.diags([dN, d1, d0, d1, dN], [-Nx, -1, 0, 1, Nx], format='lil')
    
    b = np.zeros(nunk)
    
    # 2. Aplicação do Contorno e Fonte
    for i in range(Nx):
        for j in range(Ny):
            Ic = ij2n(i, j, Nx)
            if i == 0:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TL
            elif i == Nx-1:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TR
            elif j == 0:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TB[i]
            elif j == Ny-1:
                A[Ic, :] = 0; A[Ic, Ic] = 1; b[Ic] = TT[i]
            else:
                if fonte is not None:
                    b[Ic] = fonte * h**2
                    
    A = A.tocsr()
    
    # 3. Preparação do Gauss-Seidel Esparso
    M = sparse.tril(A, format='csr')
    x = np.zeros(nunk)
    F = A.dot(x) - b
    
    k_iter = 0
    frames = [] if animation else None
    t0 = time.time()
    
    # 4. Loop Iterativo
    while np.linalg.norm(F, np.inf) > TOL and k_iter < MAXIT:
        d = spsolve_triangular(M, -F, lower=True) 
        x = x + d
        F = A.dot(x) - b
        k_iter += 1

        if animation and (k_iter % frame_skip == 0):
            frames.append(x.reshape((Ny, Nx)).copy())

    tempo = time.time() - t0
    T_grid = x.reshape((Ny, Nx))
    return T_grid, k_iter, tempo, frames

# FUNÇÃO AnimacaoTemperatura --------------------------------------------------------------------------

def AnimacaoTemperatura(frames, Nx, Ny, Lx, Ly):
    
    if frames is None or len(frames) == 0:
        print("Nenhum frame foi gerado!")
        return
    
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    X, Y = np.meshgrid(x, y)
    
    fig, ax = plt.subplots(figsize=(6,6))
    
    # 1. Encontra a temperatura mais fria e mais quente reais da simulação
    vmin_real = min(frame.min() for frame in frames)
    vmax_real = max(frame.max() for frame in frames)
    
    # 2. A CORREÇÃO: Arredonda os limites para o múltiplo de 5 mais próximo
    # Ex: se o mínimo for 12, cai para 10. Se o máximo for 41.3, sobe para 45.
    vmin = int(np.floor(vmin_real / 5) * 5)
    vmax = int(np.ceil(vmax_real / 5) * 5)
    
    # Cria os níveis de cor usando os limites bonitos e arredondados
    levels = np.linspace(vmin, vmax, 21)
    
    Z_init = frames[0].reshape(Ny, Nx)
    im = ax.contourf(X, Y, Z_init, levels=levels, cmap='jet')
    
    # 3. Configura a barra de cores
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Temperatura (°C)')
    
    # 4. Força os textos da barra a pularem exatamente de 5 em 5 graus
    ticks_exatos = np.arange(vmin, vmax + 1, 5)
    cbar.set_ticks(ticks_exatos)
    cbar.ax.set_yticklabels([f"{val}°C" for val in ticks_exatos]) # Adiciona o símbolo de grau
    
    def update(frame):
        ax.clear()
        
        Z = frame.reshape(Ny, Nx)
        
        ax.contourf(X, Y, Z, levels=levels, cmap='jet')
        ax.contour(X, Y, Z, levels=levels, linewidths=0.25, colors='k')
        
        ax.set_title(f"Evolução da Temperatura ({Nx}×{Ny})")
        ax.set_xlabel('x (m)')
        ax.set_ylabel('y (m)')
        ax.set_aspect('equal')
        
        ax.set_xticks([0, Lx/2, Lx])
        ax.set_yticks([0, Ly/2, Ly])
        
        return []
    
    ani = animation.FuncAnimation(
        fig,
        update,
        frames=frames,
        interval=200,
        repeat=False
    )
    
    plt.tight_layout() 
    plt.show()
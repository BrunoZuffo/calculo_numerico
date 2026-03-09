import numpy as np

def Assembly(conec: list[list], C: list):
    conec = np.array(conec) - 1 # estamos convertendo conec em uma matriz numpy e fazendo decrementações em cada um dos indices para ficar mais pythonic
    nv = conec.max() + 1 #pegando valor maximo dos nos e somando 1 por conta do indice reduzido anteriormente
    nc = len(conec) #simplesmente o tanto de linhas na matriz conec
    A = np.zeros(shape=(nv,nv)) #criando matriz com valores 0 em todas as posições
    for k in range(nc):
        n1 = conec[k,0] #n1 e n2 recebem os nós de cada cano, cano=k ( de 0 até numero de canos-1 )
        n2 = conec[k,1]
        
        A[n1,n1] += C[k]
        A[n1,n2] -= C[k]
        A[n2,n1] -= C[k]
        A[n2,n2] += C[k]

    return A

def SolveNetwork(conec: list[list], C:list, natm, nB, QB):

    natm = natm - 1
    nB = nB - 1
    
    Atilde = Assembly(conec,C) #recebe a matriz A
    
    Atilde[natm, :] = 0 #todas as colunas com linha natm recebem 0
    Atilde[natm, natm] = 1 

    b = np.zeros(Atilde.shape(0)) #shape retorna uma tupla com o numero de linhas e colunas, nesse caso o 0 pega o primeiro, o numero de linhas
    b[nB]=QB
    
    pressure = np.linalg.solve(Atilde, b)
    return pressure
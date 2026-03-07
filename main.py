import numpy as np

def Assembly(conec: list[list], C: list[list]):
    conec = np.array(conec) - 1 # estamos convertendo conec em uma matriz numpy e fazendo decrementações em cada um dos indices para ficar mais pythonic
    nv = conec.max() + 1 #pegando valor maximo dos nos e somando 1 por conta do indice reduzido anteriormente
    nc = len(conec) #simplesmente o tanto de linhas na matriz conec
    A = np.zeros(shape=(nv,nv)) #criando matriz com valores 0 em todas as posições
    for k in range(nc):
        n1 = conec[k,0]
        n2 = conec[k,1]
        

    return A
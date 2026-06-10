import numpy as np


def compara_arquivos_matriz(arquivo_a, arquivo_b):
    """
    Compara dois arquivos de texto contendo matrizes numéricas.

    Retorna 1 se forem iguais e -1 se forem diferentes.
    Se houver diferença, imprime onde ela apareceu primeiro.
    """
    try:
        matriz_a = np.loadtxt(arquivo_a)
        matriz_b = np.loadtxt(arquivo_b)
    except OSError as erro:
        print(f"Erro ao abrir arquivo: {erro}")
        return -1
    except ValueError as erro:
        print(f"Erro ao ler os dados numéricos: {erro}")
        return -1

    if matriz_a.shape != matriz_b.shape:
        print(f"Arquivos diferentes: formato {matriz_a.shape} != {matriz_b.shape}")
        return -1

    iguais = np.array_equal(matriz_a, matriz_b)
    if iguais:
        return 1

    diferencas = np.argwhere(matriz_a != matriz_b)
    i, j = diferencas[0]
    valor_a = matriz_a[i, j]
    valor_b = matriz_b[i, j]
    print(
        f"Primeira diferença em linha {i + 1}, coluna {j + 1}: "
        f"{arquivo_a} = {valor_a} | {arquivo_b} = {valor_b}"
    )
    return -1


if __name__ == "__main__":
    resultado_k = compara_arquivos_matriz("K.txt", "K1.txt")
    resultado_m = compara_arquivos_matriz("M.txt", "M1.txt")

    print(f"Resultado K: {resultado_k}")
    print(f"Resultado M: {resultado_m}")
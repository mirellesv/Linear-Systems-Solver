# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:30:24 2023

@author: mirelle
"""

import numpy as np

def jacobi(A, b, x0, tol, iteracoes):
    n = len(A) # n corresponde ao número de linhas
    x = np.zeros(n, dtype = 'double') # x armazena a solução atual
    xant = x0 # xant armazena a solução da iteração anterior
    for k in range(iteracoes): # Iteração do método
        for i in range(n): # Iteração feita para cada incógnita
            x[i] = b[i] # Termo constante
            for j in range(i): # Incógnitas anteriores
                x[i] -= A[i, j]*xant[j]
            for j in range(i + 1, n): # Incógnitas posteriores
                x[i] -= A[i, j]*xant[j]
            x[i] /= A[i, i] # Divisão pelo coeficiente da incógnita atual
        erro = np.linalg.norm(x - xant, np.inf) # Cálculo do erro
        print("Iteração {k:3d}: ".format(k = k+1) +
              "x = {x}, ".format(x = np.round(x, 8)) + 
              "Erro = {e:+5.8f}".format(e = erro))
        if(erro < tol): # Teste que verifica se o erro é menor que a tolerância
            return x
        xant = np.copy(x) # A solução atual é copiada ao passar a ser a anterior para a próxima iteração

# A linha pivô em M é procurada na coluna col a partir da linha lin
def linha_pivo(mat, lin, col):
    maior = abs(M[lin, col]) # O primeiro elemento é o maior valor
    
    lin_pivo = lin # Linha pivô é a linha do primeiro elemento
    
    n = len(M) # n corresponde ao número de linhas de M
    
    for i in range(lin, n): # A coluna col é percorrida a partir de lin
        if abs(M[i, col]) > maior:
            maior = abs(M[i, col]) # O pivô é atualizado
            lin_pivo = i
    return lin_pivo
    
# lin1 é trocada com lin2 em M
def troca_linha(M, lin1, lin2):
    if lin1 != lin2: # As duas linhas são diferentes?
        print('Troca de linhas: ', lin1, '<->', lin2)
        aux = np.copy(M[lin1, :])
        M[lin1, :] = np.copy(M[lin2, :])
        M[lin2, :] = aux
        print(M)
        
# A matriz M no formato de diagonal superior é resolvida
def resolve_diag_sup(M):
    n = len(M) # n corresponde ao número de linhas de M
        
    b = np.copy(M[:, n]) # b é o vetor de termos constantes
        
    x = np.arange(n, dtype = 'double') # Vetor x é criado para guardar a solução
        
    if(M[n - 1, n - 1] != 0): # Verifica se não existem zeros na diagonal principal
                              # Se não existe, o cálculo continua
        x[n - 1] = b[n - 1] / M[n - 1, n - 1] # A última linha já está isolada
        
        for i in range(n - 2, -1, -1): # As linhas são percorridas em ordem decrescente igonrando a última
            soma = 0
                
            for j in range(i + 1, n): # As incógnitas já resolvidas são somadas (depois da diagonal)
                soma += x[j] * M[i, j] # x_i = (b_i - (soma das incógnitas)) / M_i, i
                
                x[i] = (b[i] - soma)/ M[i, i]
        return x
    
    else: # Caso exista, não há soluções para o sistema
        print('Não há soluções para o sistema, dado que existem zeros na diagonal principal')
        return ""

# Eliminação gaussiana com pivotamento
def gauss_pivo(M):
    n = len(M) # n corresponde ao número de linhas de M
    
    for c in range(n - 1):
        print('\n\nColuna', c)
        l = linha_pivo(M, c, c)
        troca_linha(M, l, c)
        pivo = M[c, c]
        for l in range(c+1, n):
            print('\nL{l} <- L{l} - L{c} * '.format(l = l, c = c) +
                  '{b} / {p}'.format(b = M[l, c], p = pivo))
            M[l, :] = M[l, :] - M[c, :] * M[l, c] / pivo
            print(M)
            
        for l in range(n):
            for c in range(n):
                if(l == c) and M[l][c] == 0:
                    print('Não existem soluções para o sistema!')
                    return None

    return resolve_diag_sup(M)

def gauss_seidel(A, b, x0, tol, iteracoes):
    n = len(A)
    x = np.zeros(n, dtype = 'double')
    xant = x0
    for k in range(iteracoes):
        for i in range(n):
            x[i] = b[i]
            for j in range(i):
                x[i] -= A[i, j] * x[j]
            for j in range(i + 1, n):
                x[i] -= A[i, j]*xant[j]
            x[i] /= A[i, i]
        erro = np.linalg.norm(x - xant, np.inf)
        print('Iteração {k:3d}:'.format(k = k+1) +
              ' x = {x}, '.format(x = np.round(x, 8)) +
              'Erro = {e:5.8f}'.format(e = erro))
        if(erro < tol):
            return x
        
        xant = np.copy(x)

mat = []
M = []
n = int(input('Informe o número de variáveis do sistema: '))
print('Informe a matriz estendida (linha por linha)')
for i in range(n):
    lin = input('Linha {i}: '.format(i = i))
    lin_num = lin.split(',')
    lin_num = [float(i) for i in lin_num]
    M.append(lin_num)

b = []

for l in range(n):
    for c in range(n + 1):
        if (c == n):
            b.append(M[l][c])
            
for l in range(n):
    for c in range(n):
        mat.append(M[l][c])

M = np.array(M, dtype = 'double')
b = np.array(b, dtype = 'double')
mat = np.array(mat, dtype = 'double')

matCoef = np.array_split(mat, n)

matCoef = np.array(matCoef, dtype = 'double')

x = gauss_pivo(M)

if not x is None:
    x0 = np.array([1, 1, -1], dtype = 'double')
    x = gauss_seidel(M, b, x0, 0.0001, 20)
    print('\nSolução aproximada encontrada')
    print('x = ', x)
    print('Cálculo do raio espectral: ')
    D = np.diag(np.diag(matCoef))
    L = np.tril(matCoef) - D
    U = np.triu(matCoef) - D
    T = -np.linalg.inv(L + D).dot(U)
    C = np.linalg.inv(L + D).dot(b)
    av, _ = np.linalg.eig(T)
    raio_espectral = max(abs(av))
    print(raio_espectral)
else:
    print('Solução vazia!')

"""
x0 = np.array([1, 1, -1], dtype = 'double')

x = jacobi(M, b, x0, 0.0001, 20)

print('\nSolução aproximada encontrada: ')
print('x = ', x)"""
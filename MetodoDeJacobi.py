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

M = np.array(M, dtype = 'double')
b = np.array(b, dtype = 'double')


b = np.array(b, dtype = 'double')

print(M)
print(b)

x0 = np.array([1, 1, -1], dtype = 'double')

x = jacobi(M, b, x0, 0.0001, 20)

print('\nSolução aproximada encontrada: ')
print('x = ', x)
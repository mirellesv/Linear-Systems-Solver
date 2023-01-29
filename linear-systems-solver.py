# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 10:30:24 2023

Disciplina: Cálculo Numérico
Trabalho: Resolvedor de Sistemas Lineares
Nome: Mirelle Silva Vieira
RA: 0059636

Objetivos:
- O programa soluciona sistemas lineares com n variáveis e n equações, inseridos pelo usuário;

- O programa solicita ao usuário todas as informações do sistema (número de variáveis, matriz estendida 
  do sistema, tolerância e número máximo de iterações);

- Após a entrada dos dados, o programa mostra o sistema para que o usuário possa conferir e confirmar se 
  deseja prosseguir;

- Depois da confirmação, o resolvedor realiza um pivotamento na matriz estendida (isto evita erros quando 
  a diagonal principal possui zeros, isto é, quando o sistema não possui soluções).

- Em seguida, o critério do raio espectral é utilizado para testar se o sistema original convergirá;

- Por fim, o programa deve calcular a solução utilizando um método iterativo (Jacobi ou Gauss-Seidel),
  inserido pelo usuário;
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
        
    x[n - 1] = b[n - 1] / M[n - 1, n - 1] # A última linha já está isolada
        
    for i in range(n - 2, -1, -1): # As linhas são percorridas em ordem decrescente igonrando a última
        soma = 0
                
        for j in range(i + 1, n): # As incógnitas já resolvidas são somadas (depois da diagonal)
            soma += x[j] * M[i, j] # x_i = (b_i - (soma das incógnitas)) / M_i, i
                
            x[i] = (b[i] - soma)/ M[i, i]
    return x

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

    for l in range(n): # Se a matriz possuir algum 0 na diagonal principal,
                       # não existirá soluções para o sistema
        for c in range(n):
            if(l == c) and M[l][c] == 0:
                print('Não existem soluções para o sistema!')
                return None
            
    return resolve_diag_sup(M)

def gauss_seidel(A, b, x0, tol, iteracoes):
    n = len(A) # n é o número de linhas
    x = np.zeros(n, dtype = 'double') # Solução atual
    xant = x0 # Solução da iteração anterior
    for k in range(iteracoes): # Iterações do método
        for i in range(n): # Iterações para cada incógnita
            x[i] = b[i] # Termo constante
            for j in range(i): # Incógnitas anteriores
                x[i] -= A[i, j] * x[j]
            for j in range(i + 1, n): # Incógnitas posteriores
                x[i] -= A[i, j]*xant[j]
            x[i] /= A[i, i] # Divisão pelo coeficiente da incógnita atual
        erro = np.linalg.norm(x - xant, np.inf) # Cálculo do erro
        print('Iteração {k:3d}:'.format(k = k+1) +
              ' x = {x}, '.format(x = np.round(x, 8)) +
              'Erro = {e:5.8f}'.format(e = erro))
        if(erro < tol): # Testa se erro é menor que a tolerância
            return x
        
        xant = np.copy(x) # Solução atual é copiada para ser a anterior na próxima iteração

print('-' * 30)
print('RESOLVEDOR DE SISTEMAS LINEARES')
print('-' * 30)
while True:
    mat = [] # Matriz auxiliar para obter a matriz dos coeficientes
    M = []
    n = int(input('Informe o número de variáveis do sistema: '))
    print('Informe a matriz estendida (linha por linha)')
    for i in range(n): # Laço de repetição que permite a criação da matriz a partir
                       # dos dados inseridos pelo usuário
        lin = input('Linha {i}: '.format(i = i))
        lin_num = lin.split(',')
        lin_num = [float(i) for i in lin_num]
        M.append(lin_num)
    
    tol = float(input('Informe a tolerância: '))
    iteracoes = int(input('Informe o número máximo de iterações: '))
    
    print('-' * 30)
    
    print('Número de variáveis do sistema: ', n)
    print('Matriz estendida inserida: ')
    for l in range(n):
        print(M[l])
    print('Tolerância: ', tol)
    print('Número máximo de iterações:', iteracoes)
    resp = input('Deseja continuar (S/N)?')
    
    if(resp.lower() == 's'):
        break

b = [] # b corresponde a matriz dos termos independentes

for l in range(n):
    for c in range(n + 1):
        if (c == n): # Os termos independentes estão na coluna n
            b.append(M[l][c]) # Eles são inseridos em b
            
for l in range(n):
    for c in range(n):
        mat.append(M[l][c]) # mat será a matriz M sem os termos independentes

M = np.array(M, dtype = 'double')
b = np.array(b, dtype = 'double')
mat = np.array(mat, dtype = 'double')

matCoef = np.array_split(mat, n) # Os números são divididos em conjunto de acordo com
                                 # com as linhas onde estão na matriz

matCoef = np.array(matCoef, dtype = 'double')

x = gauss_pivo(M) # Escalonamento da matriz

print()
if not x is None: # Se não existirem zeros na diagonal principal ...
    print('Não existem zeros na diagonal principal. Logo, o sistema possui solução.')
    print('1. Método de Jacobi')
    print('2. Método de Gauss-Seidel')
    op = int(input('Selecione o método iterativo desejado: '))
    if(op == 1):
        print()
        print('Cálculo do raio espectral: ', end = '')
        D = np.diag(np.diag(matCoef))
        L = np.tril(matCoef) - D
        U = np.triu(matCoef) - D
        T = -np.linalg.inv(D).dot(L + U)
        C = np.linalg.inv(D).dot(b)
        av, _ = np.linalg.eig(T)
        raio_espectral = max(abs(av))
        print(raio_espectral)
        if(raio_espectral <= 1):
            print('Raio espectral inferior ou igual a 1. Logo, o método Jacobi converge.')
            print()
            x0_list = [] # Lista que irá armazenar os valores iniciais para todas as variáveis
            for i in range(n):
                x0_list.append(0) # Todas as variáveis utilizarão zero como valor inicial
            x0 = np.array(x0_list, dtype = 'double')
            x = jacobi(matCoef, b, x0, tol, iteracoes)
            if not x is None: # Se o número de iterações for suficiente...
                print('\nSolução aproximada encontrada')
                print('x = ', x)
            else:
                print(iteracoes, 'iterações não foram suficientes para encontrar as soluções!')
        else:
            print('Raio espectral superior a 1. Logo, o método Jacobi não converge.')
    elif(op == 2):
        print('Cálculo do raio espectral: ', end = '')
        D = np.diag(np.diag(matCoef))
        L = np.tril(matCoef) - D
        U = np.triu(matCoef) - D
        T = -np.linalg.inv(L + D).dot(U)
        C = np.linalg.inv(L + D).dot(b)
        av, _ = np.linalg.eig(T)
        raio_espectral = max(abs(av))
        print(raio_espectral)
        if(raio_espectral <= 1):
            print('Raio espectral inferior ou igual a 1. Logo, o método Gauss-Seidel converge.')
            print()
            x0_list = [] # Lista que irá armazenar os valores iniciais para todas as variáveis
            for i in range(n):
                x0_list.append(0) # Todas as variáveis utilizarão zero como valor inicial
            x0 = np.array(x0_list, dtype = 'double')
            x = gauss_seidel(matCoef, b, x0, tol, iteracoes)
            if not x is None: # Se o número de iterações for suficiente...
                print('\nSolução aproximada encontrada')
                print('x = ', x)
            else:
                print(iteracoes, 'iterações não foram suficientes para encontrar as soluções!')
        else:
            print('Raio espectral superior a 1. Logo, o método Gauss-Seidel não converge.')
else:
    print('Solução vazia!')
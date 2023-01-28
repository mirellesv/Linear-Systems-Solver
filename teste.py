# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 14:20:55 2023

@author: mirelle
"""

import numpy as np

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
    for c in range(n):
        mat.append(M[l][c])

for l in range(n):
    for c in range(n + 1):
        if (c == n):
            b.append(M[l][c])

M = np.array(M, dtype = 'double')
b = np.array(b, dtype = 'double')
mat = np.array(mat, dtype = 'double')

matCoef = np.array_split(mat, 3)

matCoef = np.array(matCoef, dtype = 'double')

print(matCoef)


D = np.diag(np.diag(matCoef))
L = np.tril(matCoef) - D
U = np.triu(matCoef) - D
TG = -np.linalg.inv(L + D).dot(U)
CG = np.linalg.inv(L + D).dot(b)
av,_ = np.linalg.eig(TG)
raio_espectral = max(abs(av))
print(raio_espectral)
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 28 20:01:27 2023

@author: mirelle
"""

import numpy as np

a = np.array([[3, 1, -1], [-1, -4, 1], [1, -2, -5]], dtype = 'double')
b = np.array([2, -10, 10], dtype = 'double')
D = np.diag(np.diag(a))
L = np.tril(a) - D
U = np.triu(a) - D
TG = -np.linalg.inv(L + D).dot(U)
CG = np.linalg.inv(L + D).dot(b)
av,_ = np.linalg.eig(TG)
raio_espectral = max(abs(av))
print(raio_espectral)
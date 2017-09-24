# -*- coding: utf-8 -*-
"""
Created on Sat Sep 23 15:04:10 2017

@author: roco33
"""



import numpy as np
import math

S0 = 100
K = 105
r = 0.01
div = 0
sigma = 0.3
T = 1
N = 5000

delta = T/N
u = math.exp(sigma * math.sqrt(delta))
d = 1/u
qu = (math.exp(r*delta) - d)/(u - d)
qd = 1 - qu

VStock = np.zeros((25000,25000))
VOption = np.zeros((25000,25000))

j = N
for i in range(j+1):
    VStock[j,i] = S0 * u**i * d ** (N - i)
    VOption[j,i] = max(K - VStock[j,i], 0)

for j in range(N-1,-1,-1):
    for i in range(j,-1,-1):
        VOption[j,i] = max(math.exp(-r * delta) * (qu * VOption[j+1,i+1] + qd * 
               VOption[j+1,i]),0)

print(VOption[0,0])

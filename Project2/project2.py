# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 03:05:48 2017

@author: roco3
"""

import numpy as np

jmax = 100
imax = 1250
s0 = 156.99
D = 15.875
T = 1.25
r = 0.02
sig = 0.2
div = 0.02
smin = 0
smax = 330
delt = T / imax
dels = (smax - smin) / jmax


s = np.zeros(jmax + 1)
A = np.zeros(jmax + 1)
B = np.zeros(jmax + 1)
C = np.zeros(jmax + 1)
D = np.zeros(jmax + 1)
V = np.zeros(jmax + 1)


# triggered
VT = np.zeros(jmax + 1)
for j in range(jmax + 1):
    s[j] = smin + j * dels
    
    if s[j] / s0 < 0.7:
        VT[j] = s[j] / s0 * 1000
    elif (s[j] / s0 >= 0.7) & (s[j] / s0 < 1):
        VT[j] = s[j] / s0 * 1000 + D
    else:
        VT[j] = 1000 + D
    # A, B, C [0 to jmax]
    A[j] = 0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j
    B[j] = -1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ** 2
    C[j] = 0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j

for i in range(imax - 1, -1, -1):
    # lower boundary
    A[0] = 0
    B[0] = 1
    C[0] = 0
    D[0] = B[0] * V[0] + C[0] * VT[1] 
    # D across 1 to (jmax-1)
    for j in range(1, jmax, 1):
        if i % (imax / 5) == 0:
            if s[j] >= s0:
                VT[j] = 1000
            elif s[j] / s0 >= 0.7 and s[j] / s0 < 1:
                VT[j] = VT[j] + D
        D1 = -VT[j - 1] * A[j]
        D2 = -VT[j] * B[j]
        D3 = -VT[j + 1] * C[j]
        D[j] = D1 + D2 + D3
    A[jmax] = 0
    B[jmax] = 1
    C[jmax] = 0
    n = int(i / (imax / 5)) * (imax / 5) # the last review period 
    D[jmax] = (1000 + D) * np.exp(-div * ((i - n) * delt))
    
    VT = LU(jmax, A, B, C, D, VT)

# not triggered
VNT = np.zeros(jmax + 1)
for j in range(70, jmax + 1, 1):
    s[j] = smin + j * dels 
    
    if s[j] / s0 < 0.7:
        VNT[j] = 1000
    else:
        VNT[j] = 1000 + D
    A[j] = 0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j
    B[j] = -1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ** 2
    C[j] = 0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j

for i in range(jmax + 1, -1, -1):
    # lower boundary
    A[0] = 0
    B[0] = 1
    C[0] = 0
    D[0] = B[0] * V[0] + C[0] * VNT[1]
    # D accross 1 to (jmax+1)
    for j in range(1, jmax, 1):
        if i % (imax / 5) == 0:
            if s[j] >= s0:
                VNT[j] = 1000
            elif s[j] / s0 >= 0.7 and s[j] / s0 < 1:
                VNT[j] = VNT[j] + D
        D1 = VNT[j - 1] * A[j]
        D2 = VNT[j] * B[j]
        D3 = VNT[j + 1] * C[j]
        D[j] = D1 + D2 + D3
    A[jmax] = 0
    B[jmax] = 1
    C[jmax] = 0
    n = int(i / (imax / 5)) * (imax / 5) # the last review period 
    D[jmax] = (1000 + D) * np.exp(-div * ((i - n) * delt))
    
    VNT = LU(jmax, A, B, C, D, VNT)






def LU(jmax, A, B, C, D, V):
    alpha = np.zeros(jmax + 1)
    S = np.zeros(jmax + 1)
    alpha[0] = B[0]
    S[0] = D[0]
    for j in range(1,jmax + 1, 1):
        alpha[j] = B[j] - A[j] * C[j - 1] /alpha[j - 1]
        S[j] = D[j] - A[j] / alpha[j - 1] * S[j - 1]
    V[jmax] = S[jmax] / alpha[jmax]
    for j in range(jmax - 1, -1, -1):
        V[j] = 1 / alpha[j] * (S[j] - C[j] * V[j + 1])
    return V
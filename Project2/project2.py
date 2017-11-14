# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 03:05:48 2017

@author: roco3
"""

import numpy as np

imax = 1250
s0 = 156.99
smin = 0
smax = 3 * s0
TAC = 100 #auto call index, should be multiplier of 10
dels = s0 / TAC # delta s
TR = 0.7 * s0 / dels # trigger index
jmax = 3 * TAC 
cou = 15.875 # coupon
T = 1.25 # maturity
r = 0.02 
sig = 0.2
div = 0.02

delt = T / imax


s = np.zeros(jmax + 1)
A = np.zeros(jmax + 1)
B = np.zeros(jmax + 1)
C = np.zeros(jmax + 1)
D = np.zeros(jmax + 1)
V = np.zeros((imax + 1, jmax + 1))
VT = np.zeros((imax + 1, jmax + 1))


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


# triggered

# at maturity
for j in range(jmax + 1):
    s[j] = smin + j * dels
    
    if s[j] / s0 < 0.7:
        VT[imax+1,j] = s[j] / s0 * 1000
    elif (s[j] / s0 >= 0.7) and (s[j] / s0 < 1):
        VT[imax+1,j] = s[j] / s0 * 1000 + cou
    else:
        VT[imax+1,j] = 1000 + cou
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
        D1 = -VT[i+1,j - 1] * A[j]
        D2 = -VT[i+1,j] * B[j]
        D3 = -VT[i+1,j + 1] * C[j]
        D[j] = D1 + D2 + D3
        if i % (imax / 5) != 0: # not review date
            A[jmax] = 0    
            B[jmax] = 1
            C[jmax] = 0
            n = int(i / (imax / 5)) * (imax / 5) # the last review period 
            D[jmax] = (1000 + cou) * np.exp(-(r - div) * ((i - n) * delt)) 
            VT[i,:] = LU(jmax, A, B, C, D, VT)
        elif i % (imax / 5) == 0: # review date
            VT[TAC] = 1000
            A[TAC] = 0
            B[TAC] = 1
            C[TAC] = 0
            D[TAC] = 1000
            VT[0: TAC + 1] = LU(TAC, A[0: TAC + 1], B[0: TAC + 1], C[0: TAC + 1], D[0: TAC + 1], VT[0: TAC + 1])
            for j in range(TAC + 1, jmax + 1, 1):
                VT[j] = 1000          
                    

# not triggered
V = np.zeros(jmax + 1)
# TR to jmax
for j in range(TR, jmax + 1, 1):
    s[j] = smin + j * dels 
    
    if s[j] / s0 < 0.7:
        V[j] = 1000
    else:
        V[j] = 1000 + cou
    
for i in range(imax - 1, -1, -1):
    # lower boundary
    A[0] = 0
    B[0] = 1
    C[0] = 0
    D[0] = B[0] * V[0] + C[0] * V[1]
    # D accross 1 to (jmax+1)
    for j in range(1, jmax, 1):
        if i % (imax / 5) == 0:
            if s[j] >= s0:
                V[j] = 1000
            elif s[j] / s0 >= 0.7 and s[j] / s0 < 1:
                V[j] = V[j] + cou
        D1 = V[j - 1] * A[j]
        D2 = V[j] * B[j]
        D3 = V[j + 1] * C[j]
        D[j] = D1 + D2 + D3
    A[jmax] = 0
    B[jmax] = 1
    C[jmax] = 0
    n = int(i / (imax / 5)) * (imax / 5) # the last review period 
    D[jmax] = (1000 + cou) * np.exp(-(r-div) * ((i - n) * delt))
    
    VNT = LU(jmax, A, B, C, D, V)







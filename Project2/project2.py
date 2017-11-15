# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 03:05:48 2017

@author: roco3
"""

import numpy as np
import matplotlib.pyplot as plt

imax = 1250 # total period
RD = [0, 0.25, 0.5, 0.75, 1, 1.25]
RD = [i * imax / 1.25 for i in RD] # review date
s0 = 156.99
smin = 0
smax = 3 * s0
TAC = 100 #auto call index, should be multiplier of 10
dels = s0 / TAC # delta s
TR = round(0.7 * s0 / dels) # trigger index
jmax = 3 * TAC 
cou = 15.875 # coupon
T = 1.25 # maturity
r = 0.02 
sig = 0.18
div = 0.022

delt = T / imax


A = np.zeros(jmax + 1)
B = np.zeros(jmax + 1)
C = np.zeros(jmax + 1)
D = np.zeros(jmax + 1)
V = np.zeros((imax, jmax + 1))
VT = np.zeros((imax, jmax + 1))

s = np.array([smin + j * dels for j in range(jmax + 1)]) # stock price grid


def LU(A, B, C, D, X):
    jmax = len(A)-1
    alpha = np.zeros(jmax + 1)
    Y = np.zeros(jmax + 1)
    S = np.zeros(jmax + 1)
    alpha[0] = B[0]
    S[0] = D[0]
    for j in range(1,jmax + 1, 1):
        alpha[j] = B[j] - A[j] * C[j - 1] /alpha[j - 1]
        S[j] = D[j] - A[j] / alpha[j - 1] * S[j - 1]
    Y[jmax] = S[jmax] / alpha[jmax]
    for j in range(jmax - 1, -1, -1):
        Y[j] = 1 / alpha[j] * (S[j] - C[j] * X[j + 1])
    return Y


# triggered

# at maturity
for j in range(jmax + 1):
    
    if s[j] / s0 < 0.7:
        VT[imax -1 ,j] = s[j] / s0 * 1000
    elif (s[j] / s0 >= 0.7) and (s[j] / s0 < 1):
        VT[imax -1 ,j] = s[j] / s0 * 1000 + cou
    else:
        VT[imax -1 ,j] = 1000 + cou


for i in range(imax - 1, 0, -1):
    for j in range(jmax+1):
        # A, B, C [0 to jmax]
        A[j] = 0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j
        B[j] = -1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ** 2
        C[j] = 0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j
    # lower boundary
    A[0] = 0
    B[0] = 1
    C[0] = 0
    D[0] = 0
    # D across 1 to (jmax-1)
    for j in range(1, jmax, 1):
        D1 = VT[i,j - 1] * A[j]
        D2 = VT[i,j] * B[j]
        D3 = VT[i,j + 1] * C[j]
        D[j] = D1 + D2 + D3
        
    if i not in RD: # regular day
        A[jmax] = 0    
        B[jmax] = 1
        C[jmax] = 0
        n = [t for t in RD if t <= i][-1] # the last review period 
        D[jmax] = (1000 + cou) * np.exp(-(r - div) * (i - n) * delt) 
        VT[i-1,:] = LU(A, B, C, D, VT[i,:])
        
    elif i in RD: # review date
        VT[i-1,TAC] = 1000
        A[TAC] = 0
        B[TAC] = 1
        C[TAC] = 0
        D[TAC] = 1000
        VT[i-1,0: TAC] = LU(A[0: TAC], B[0: TAC], C[0: TAC], D[0: TAC], VT[i, 0: TAC])
        for j in range(TAC + 1, jmax + 1, 1):
            VT[i-1,j] = 1000
        for j in range(TR, jmax + 1, 1):
            VT[i-1,j] = VT[i-1,j] + cou
                    

A = np.zeros(jmax + 1)
B = np.zeros(jmax + 1)
C = np.zeros(jmax + 1)
D = np.zeros(jmax + 1)
# not triggered

# TR to jmax
for j in range(TR, jmax + 1, 1):
    
    if s[j] / s0 < 0.7:
        V[imax-1,j] = 1000
    else:
        V[imax-1,j] = 1000 + cou

    
for i in range(imax - 1, 0, -1):
    for j in range(TR, jmax + 1, 1):
        A[j] = 0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j
        B[j] = -1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ** 2
        C[j] = 0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j
    # lower boundary
    V[i-1,TR] = VT[i-1,TR]
    A[TR] = 0
    B[TR] = 1
    C[TR] = 0
    D[TR] = V[i-1,j]
    # D accross 1 to (jmax+1)
    for j in range(TR+1, jmax, 1):
        D1 = V[i,j - 1] * A[j]
        D2 = V[i,j] * B[j]
        D3 = V[i,j + 1] * C[j]
        D[j] = D1 + D2 + D3
		
    if i not in RD: # regular day
        A[jmax] = 0
        B[jmax] = 1
        C[jmax] = 0
        n = [t for t in RD if t <= i][-1] # the last review period 
        D[jmax] = (1000 + cou) * np.exp(-(r - div) * (i - n) * delt) 
        V[i-1,TR:jmax] = LU(A[TR:jmax], B[TR:jmax], C[TR:jmax], D[TR:jmax], V[i,TR:jmax])
	
    elif i in RD: # review day
        V[i-1,TAC] = 1000
        A[TAC] = 0
        B[TAC] = 1
        C[TAC] = 0
        D[TAC] = 1000
        V[i-1,TR:TAC] = LU(A[TR:TAC], B[TR:TAC], C[TR:TAC], D[TR:TAC], V[i, TR:TAC])
        for j in range(TAC + 1, jmax + 1,1):
            V[i-1,j] = 1000
        for j in range(TR, jmax + 1, 1):
            V[i-1,j] = V[i-1,j] + cou
	




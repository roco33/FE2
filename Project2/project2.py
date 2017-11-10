# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 21:09:55 2017

@author: roco33
"""

'''
Value European call using finite difference method
test error through different jmax
'''

import numpy as np

Input = {
        'nmin' : 40, # jmax min 
        'nmax' : 500, # jmax max
        'jmaxtep' : 20, # test jmax in jmaxtep steps
        'imax' : 1000, # imax
        'S0' : 40, # stock price
        'X' : 40, # strike 
        'sig' : 0.2, # stock volatility
        'r' : 0.06, # risk-free
        'div' : 0, # dividend yield
        'T' : 0.5, # time horizon
        'smin' : 0, # lower boundary
        'smax' : 100, # upper boundary
        }


delt = T / imax # dt
dels = (smax - smin) / jmax # ds




def payoff(jmax, dels, delt, smin, sig, r, div,):
    V = np.zeros(jmax)
    A = np.zeros(jmax)
    B = np.zeros(jmax)
    C = np.zeros(jmax)
    for j in range(jmax):
        V[j] = (smin + j * dels) - X
        if V[j] < 0:
            V[j] = 0
        A[j] = 0.25 * sig ^ 2 * j ^ 2 - 0.25 * (r - div) * j
        B[j] = -1 / delt - 0.5 * r - 0.5 * sig ^ 2 * j ^ 2
        C[j] = 0.25 * sig ^ 2 * j ^ 2 + 0.25 * (r - div) * j
    return V, A, B, C


def bounaryCond():
    for i in range(1,imax+1,1):
        # lower boundary
        A[0] = 0
        B[0] = 1
        C[0] = 0
        D[0] = 0
        # D accross j levels
        for j in range(1,jmax,1):
            D1 = -V[j - 1] * (0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j)
            D2 = -V[j] * (1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ^ 2)
            D3 = -V[j + 1] * (0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j)
            D[j] = D1 + D2 + D3
        # upper boundary equations
        A[jmax] = 0
        B[jmax] = 1
        C[jmax] = 0
        D[jmax] = smax * np.exp(-div * (i * delt)) - X * np.exp(-r * (1 * delt))
    return A, B, C, D


def LU():
    alpha[0] = B[0]
    Sum[0] = D[0]
    for j in range(1, jmax + 1, 1):
        alpha[j] = B[j] - ((A[j] * C[j - 1]) / alpha[j - 1])
        Sum[j] = D[j] - ((A[j] * Sum[j - 1]) / alpha[j - 1])
    V[jmax] = Sum[jmax] / alpha[jmax]
    for j in range(jmax-1, -1, -1):
        V[j] = (Sum[j] - C[j] * V[j + 1]) / alpha[j]
    return V
        
        
def main(I):
    
    #payoff
    #boundary condition
    #   LU decomposition
    #interplolation
    return

    
            
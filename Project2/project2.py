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

nmin = 40 # jmax min 
nmax = 500 # jmax max
jmaxtep = 20 # test jmax in jmaxtep steps
imax = 1000
S0 = 40 # stock price
X = 40 # strike 
sig = 0.2 # stock volatility
r = 0.06 # risk-free
div = 0 # dividend yield
T = 0.5 # time horizon
smin = 0 # lower boundary
smax = 100 # upper boundary

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
    for i in range(imax-1,-1,-1):
        pass
    return


def main():
    #payoff
    #boundary condition
    #   LU decomposition
    #interplolation
    return

    
            
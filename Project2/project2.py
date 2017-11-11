# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 03:05:48 2017

@author: roco3
"""

import numpy as np

def LU(jmax):
    alpha = np.zeros(jmax + 1)
    S = np.zeros(jmax + 1)
    A = np.zeros(jmax + 1)
    B = np.zeros(jmax + 1)
    C = np.zeros(jmax + 1)
    D = np.zeros(jmax + 1)
    V = np.zeros(jmax + 1)
    alpha[0] = B[0]
    S[0] = D[0]
    for j in range(1,jmax + 1, 1):
        alpha[j] = B[j] - A[j] * C[j - 1] /alpha[j - 1]
        S[j] = D[j] - A[j] / alpha[j - 1] * S[j - 1]
    V[jmax] = S[jmax] / alpha[jmax]
    for j in range(jmax - 1, -1, -1):
        V[j] = 1 / alpha[j] * (S[j] - C[j] * V[j + 1])
    return V
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 11 03:05:48 2017

@author: roco3
"""

import numpy as np
import matplotlib.pyplot as plt


def LU(a, b, c, d):
    jjmax = len(a)-1
    alpha = np.zeros(jjmax + 1)
    Y = np.zeros(jjmax + 1)
    S = np.zeros(jjmax + 1)
    alpha[0] = b[0]
    S[0] = d[0]
    for j in range(1,jjmax + 1, 1):
        alpha[j] = b[j] - a[j] * c[j - 1] /alpha[j - 1]
        S[j] = d[j] - a[j] / alpha[j - 1] * S[j - 1]
    Y[jjmax] = S[jjmax] / alpha[jjmax]
    for j in range(jjmax - 1, -1, -1):
        Y[j] = 1 / alpha[j] * (S[j] - c[j] * Y[j + 1])
    return Y


def CN(TAC):
    imax = 1250 # total period
    RD = [0.25, 0.5, 0.75, 1, 1.25]
    RD = [i * imax / 1.25 for i in RD] # review date
    s0 = 156.99
    smin = 0
    # smax = 3 * s0
    # TAC = 100 #auto call index, should be multiplier of 10
    dels = s0 / TAC # delta s
    TR = round(0.7 * s0 / dels) # trigger index
    jmax = 3 * TAC 
    cou = 15.875 # coupon
    T = 1.25 # maturity
    r = 0.02 
    sig = 0.24
    div = 0.0065
    
    delt = T / imax
    
    A = np.zeros(jmax + 1)
    B = np.zeros(jmax + 1)
    C = np.zeros(jmax + 1)
    D = np.zeros(jmax + 1)
    V = np.zeros((imax, jmax + 1))
    VT = np.zeros((imax, jmax + 1))
    
    s = np.array([smin + j * dels for j in range(jmax + 1)]) # stock price grid
    
    # triggered
    
    # at maturity
    for j in range(jmax + 1):
        if s[j] / s0 < 1:
            VT[imax -1 ,j] = s[j] / s0 * 1000
        else:
            VT[imax -1 ,j] = 1000 + cou
        if s[j] / s0 >= 0.7:
            VT[imax -1 ,j] = VT[imax - 1, j] + cou
    
    
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
            D1 = -VT[i,j - 1] * (0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j)
            D2 = -VT[i,j] * (1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ** 2)
            D3 = -VT[i,j + 1] * (0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j)
            D[j] = D1 + D2 + D3
            
        if i not in RD: # regular day
            A[jmax] = 0    
            B[jmax] = 1
            C[jmax] = 0
            n = [t-1 for t in RD if t-1 >= i][0] # the last review period 
            D[jmax] = (1000 + cou) * np.exp(-(r - div) * (n - i) * delt) 
            VT[i-1,:] = LU(A, B, C, D)
        
        elif i in RD: # review date
            VT[i-1,TAC] = 1000
            A[TAC] = 0
            B[TAC] = 1
            C[TAC] = 0
            D[TAC] = 1000
            VT[i-1,0: TAC] = LU(A[0: TAC], B[0: TAC], C[0: TAC], D[0: TAC])
            for j in range(TAC, jmax + 1, 1):
                VT[i-1,j] = 1000
            for j in range(TR, jmax + 1, 1):
                VT[i-1,j] = VT[i-1,j] + cou
        else:
            print("error 1")
                        
    
    A = np.zeros(jmax + 1)
    B = np.zeros(jmax + 1)
    C = np.zeros(jmax + 1)
    D = np.zeros(jmax + 1)
    # not triggered
    
    # TR to jmax
    for j in range(TR, jmax + 1, 1):
        
        if s[j] / s0 >= 0.7:
            V[imax-1,j] = 1000 + cou
#        else:
#            V[imax-1,j] = 1000
    
        
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
        D[TR] = V[i-1,TR]
        # D accross 1 to (jmax+1)
        for j in range(TR+1, jmax, 1):
            D1 = -V[i,j - 1] * (0.25 * sig ** 2 * j ** 2 - 0.25 * (r - div) * j)
            D2 = -V[i,j] * (1 / delt - 0.5 * r - 0.5 * sig ** 2 * j ** 2)
            D3 = -V[i,j + 1] * (0.25 * sig ** 2 * j ** 2 + 0.25 * (r - div) * j)
            D[j] = D1 + D2 + D3
    		
        if i not in RD: # regular day
            A[jmax] = 0
            B[jmax] = 1
            C[jmax] = 0
            n = [t-1 for t in RD if t-1 >= i][0] # the next review period 
            D[jmax] = (1000 + cou) * np.exp(-(r - div) * (n - i) * delt) 
            V[i-1,TR:jmax] = LU(A[TR:jmax], B[TR:jmax], C[TR:jmax], D[TR:jmax])
    	
        elif i in RD: # review day
            V[i-1,TAC] = 1000
            A[TAC] = 0
            B[TAC] = 1
            C[TAC] = 0
            D[TAC] = 1000
            V[i-1,TR:TAC] = LU(A[TR:TAC], B[TR:TAC], C[TR:TAC], D[TR:TAC])
            for j in range(TAC, jmax + 1,1):
                V[i-1,j] = 1000
            for j in range(TR, jmax + 1, 1):
                V[i-1,j] = V[i-1,j] + cou
        else:
            print("error 2")

    return 3*TAC, V[0,TAC]


def main():
    Tac = [100 + k * 10 for k in range(10)]
    J = []
    VV = []
    k = 0
    for i in Tac:
        value = CN(i)
        J.append(value[0])
        VV.append(value[1])
        k = k + 1
        print(str(int(k/len(Tac)*100)) + '%')
    # plt.show()
    plt.subplot(121)
    plt.plot(J,VV)
    slope = np.array([abs(VV[i+1]-VV[i])/30 for i in range(len(VV)-1)])
    plt.subplot(122)
    plt.plot(J[:-1],np.log(slope))
    plt.show()


if __name__ == '__main__':
    main()
        
# -*- coding: utf-8 -*-

import numpy as np
import math
from scipy.stats import norm
import matplotlib.pyplot as plt


# 1.
# a. calculate the B-S value of the put

S0 = 100
K = 105
r = 0.01
div = 0
sigma = 0.3
T = 1

d1 = 1 / (math.sqrt(T)) * (math.log(S0/K) + ((r - div) + sigma**2 / 2) * T)
d2 = d1 - sigma * math.sqrt(T)

P = norm.cdf(-d2) * K * math.exp(-r*T) - norm.cdf(-d1) * S0



# binomial valuation function
def binomial(S0, K, r, div, sigma, T, N, method):
    delta = T/N
    VStock = np.zeros((601,601))
    VOption = np.zeros((601,601))
    if method == 'CCR':
        u = math.exp(sigma * math.sqrt(delta))
        d = 1/u
    elif method == 'RB':
        u = math.exp((r-div-0.5*sigma**2)*delta+sigma*math.sqrt(delta))
        d = math.exp((r-div-0.5*sigma**2)*delta-sigma*math.sqrt(delta))
    elif method == 'LR':
        d1 = 1 / (math.sqrt(T)) * (math.log(S0/K) + ((r - div) + sigma**2 / 2) * T)
        d2 = d1 - sigma * math.sqrt(T)
        k = 1 if d2 > 0 else -1
        l = 1 if d1 > 0 else -1
        q = 1/2 + k * math.sqrt(1/4-1/4*math.exp(-(d2/(N+1/3))**2*(N+1/6)))
        q_star = 1/2 + l * math.sqrt(1/4-1/4*math.exp(-(d1/(N+1/3))**2*(N+1/6)))
        u = math.exp((r-div)*delta)*q_star/q
        d = (math.exp((r-div)*delta)-q*u)/(1-q)
    else:
        print("Unknown method")
        return
    qu = (math.exp(r*delta) - d)/(u - d)
    qd = 1 - qu
    j = N
    for i in range(j+1):
        VStock[j,i] = S0 * u**i * d ** (N - i)
        VOption[j,i] = max(K - VStock[j,i], 0)
    for j in range(N-1,-1,-1):
        for i in range(j,-1,-1):
            VOption[j,i] = math.exp(-r * delta) * (qu * VOption[j+1,i+1] +
                   qd * VOption[j+1,i])

    return VOption[0,0]

# b.
# CCR

P1 = np.array([])

for N in range(50,501):
    P1=np.append(P1,binomial(S0,K,r,div,sigma,T,N,'CCR'))

plt.plot(np.arange(50,501), P1-P)
plt.show()


# c.
# R&B

P2 = np.array([])

for N in range(50,501):
    P2 = np.append(P2,binomial(S0,K,r,div,sigma,T,N,'RB'))

plt.plot(np.arange(50,501), P2-P)
plt.show()


# d.
# L&R

P3 = np.array([])

for N in range(50,501):
    P3 = np.append(P3,binomial(S0,K,r,div,sigma,T,N,'LR'))

plt.plot(np.arange(50,501), P3-P)
plt.show()


# e.
N = np.array([25, 50, 100, 150, 200, 250])

P4 = np.array([])
P5 = np.array([])

for n in N:
    P4 = np.append(P4,binomial(S0,K,r,div,sigma,T,n,'CCR'))
    P5 = np.append(P5,binomial(S0,K,r,div,sigma,T,2*n,'CCR'))

plt.plot(N,P4-P,label="simple lattic")
plt.plot(N,(P4*N-P5*2*N)/(-N)-P, label="extapolating")
plt.show()


N1 = np.array([51, 101, 151, 201, 251])

P6 = np.array([])
P7 = np.array([])

for n in N1:
    P6 = np.append(P6,binomial(S0,K,r,div,sigma,T,n,'LR'))
    P7 = np.append(P7,binomial(S0,K,r,div,sigma,T,2*(n-1)-1,'LR'))

plt.plot(N1,P6-P,N1,P7-P)
plt.show()


# 2.
# a. memory error
def CCRA(S0,K,r,div,sigma,T,N):
    delta = T/N
    u = math.exp(sigma * math.sqrt(delta))
    d = 1/u
    qu = (math.exp(r*delta) - d)/(u - d)
    qd = 1 - qu
    VStock = np.zeros((30000,30000))
    VOption = np.zeros((30000,30000))
    j = N
    for i in range(j+1):
        VStock[j,i] = S0 * u**i * d ** (N - i)
        VOption[j,i] = max(K - VStock[j,i], 0)
    for j in range(N-1,-1,-1):
        for i in range(j,-1,-1):
            VStock[j,i] = S0 * u**i *d ** (j - i)
            VOption[j,i] = max(math.exp(-r * delta) * (qu * VOption[j+1,i+1] + 
                   qd * VOption[j+1,i]),max(K - VStock[j,i],0))
    return VOption[0,0]

# initial position

n = 10000
m = 0
N = np.zeros((30000))
N[0] = n
Er = np.zeros((30000))
Er[0] = abs(CCRA(S0,K,r,div,sigma,T,n) - CCRA(S0,K,r,div,sigma,T,n-1))


    
# b

#PA1 = np.array([])
#
#for N in range(50, 501):
#    PA1.append(CCRA(S0,K,r,div,sigma,T,N))
#
#plt.plot(np.arange(50,501), PA1)
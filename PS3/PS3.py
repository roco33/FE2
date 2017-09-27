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

d1 = 1 / (sigma * math.sqrt(T)) * (math.log(S0/K) + ((r - div) + sigma**2
                            / 2) * T)
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

S0 = 100
K = 105
r = 0.01
div = 0
sigma = 0.3
T = 1


# a. memory error
def binomialA(S0,K,r,div,sigma,T,N,method):
    delta = T/N
    if method in ['CCR','BD']:
        u = math.exp(sigma * math.sqrt(delta))
        d = 1/u
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
    VStock = np.zeros((11000,11000))
    VOption = np.zeros((11000,11000))
    j = N
    if method in ['CCR','LR']:
        for i in range(j+1):
            VStock[j,i] = S0 * u**i * d ** (N - i)
            VOption[j,i] = max(K - VStock[j,i], 0)
        for j in range(N-1,-1,-1):
            for i in range(j,-1,-1):            
                VStock[j,i] = S0 * u**i *d ** (j - i)
                VOption[j,i] = max(math.exp(-r * delta) * (qu * VOption[j+1,
                       i+1] + qd * VOption[j+1,i]),max(K - VStock[j,i],0))
    elif method == 'BD':
        for i in range(j):
            VStock[j-1,i] = S0 * u**i * d**(j-1-i)
            d1 = 1 / (sigma * math.sqrt(T-delta)) * (math.log(VStock[j-1,i]
            /K) + ((r - div) + sigma**2/ 2) * (T-delta))
            d2 = d1 - sigma * math.sqrt(T-delta)
            VOption[j-1,i] = max(K - VStock[j-1,i], norm.cdf(-d2) * K * math.exp
                   (-r*(T-delta)) - norm.cdf(-d1) * VStock[j-1,i])
        for j in range(N-2,-1,-1):
            for i in range(j,-1,-1):
                VStock[j,i] = S0 * u**i * d**(j-i)
                VOption[j,i] = max(K-VStock[j,i], math.exp(-r*delta)*(qu*
                       VOption[j+1,i+1])+qd*VOption[j+1,i])
    return VOption[0,0]

# initial position

n = 10000
m = 0
PA = binomialA(S0,K,r,div,sigma,T,n,'CCR')
Er = abs(binomialA(S0,K,r,div,sigma,T,n,'CCR') - binomialA(S0,K,r,div,sigma,T,
         n-1,'CCR'))

# b

# CCR
PA1 = np.array([])

for N in range(50, 501):
    PA1= np.append(PA1, binomialA(S0,K,r,div,sigma,T,N,'CCR'))

plt.plot(np.arange(50,501), PA1-PA)
plt.show()

# B&D
PA2 = np.array([])

for N in range(50,501):
    PA2 = np.append(PA2, binomialA(S0,K,r,div,sigma,T,N,'BD'))

plt.plot(np.arange(50,501), PA2-PA)
plt.show()


# L&R
PA3 = np.array([])

for N in range(50,501):
    PA3 = np.append(PA3, binomialA(S0,K,r,div,sigma,T,N,'LR'))

plt.plot(np.arange(50,501), PA3-PA)
plt.show()



# 3

S0 = 100
K = 100
r = 0.1
div = 0
sigma = 0.3
T = 0.5
B = 95

d1 = 1 / (sigma * math.sqrt( T)) * (math.log(S0/K) + ((r - div) + sigma**2 / 
                            2) * T)
d2 = d1 - sigma * math.sqrt(T)
h1 = (math.log(B**2/(K*S0))+(r-div+1/2*sigma**2)*T)/(sigma*math.sqrt(T))
h2 = (math.log(B**2/(K*S0))+(r-div-1/2*sigma**2)*T)/(sigma*math.sqrt(T))

CB = S0 * math.exp(-div*T) * norm.cdf(d1) - K * math.exp(-r*T)*norm.cdf(d2) -\
(B/S0)**(1+2*r*sigma**(-2)) *S0* norm.cdf(h1) + (B/S0)**(-1+2*r*sigma**(-2)) *\
K * math.exp(-r*T)*norm.cdf(h2)

# a

def binomialB(S0,K,r,div,sigma,T,N,B,method):
    delta = T/N
    if method == 'CCR':
        u = math.exp(sigma * math.sqrt(delta))
        d = 1/u
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
    VStock = np.zeros((600,600))
    VOption = np.zeros((600,600))
    j = N
    for i in range(j+1):
        VStock[j,i] = S0 * u**i * d ** (N - i)
        VOption[j,i] = max(VStock[j,i] - K, 0)
    for j in range(N-1,-1,-1):
        for i in range(j,-1,-1):
            VStock[j,i] = S0 * u**i *d ** (j - i)
            VOption[j,i] = (0 if VStock[j,i] < B else math.exp(-r * delta) * 
                   (qu * VOption[j+1,i+1] + qd * VOption[j+1,i]))
    return VOption[0,0]


CB1 = np.array([])

for N in range(50,501):
    CB1 = np.append(CB1, binomialB(S0,K,r,div,sigma,T,N,B,'CCR'))

plt.plot(np.arange(50,501), CB1-CB)
plt.show()




#4
S0 = 100
K = 100
B = 95
r = 0.1
div = 0
sigma = 0.3
T = 0.2


def CRRDB(S0,K,r,div,sigma,T,N,B):
    delta = T/N
    u = math.exp(sigma * math.sqrt(delta))
    d = 1/u
    qu = (math.exp(r*delta) - d)/(u - d)
    qd = 1 - qu
    VStock = np.zeros((2000,2000))
    VOption = np.zeros((2000,2000))
    j = N
    for i in range(j+1):
        VStock[j,i] = S0 * u**i * d ** (N - i)
        VOption[j,i] = max(VStock[j,i] - K, 0)
    for j in range(N-1,-1,-1):
        for i in range(j,-1,-1):
            VStock[j,i] = S0 * u**i *d ** (j - i)
            VOption[j,i] = (0 if (VStock[j,i] < B)&((j*delta) in [0.04,0.08,
                   0.12,0.16]) else math.exp(-r * delta) * (qu * VOption[j+1,
                            i+1] + qd * VOption[j+1,i]))
    return VOption[0,0]

DCB = 5.6711051343

DCB1 = np.array([])

for N in range(50,1010,10):
    DCB1 = np.append(DCB1,CRRDB(S0,K,r,div,sigma,T,N,B))

plt.plot(np.arange(50,1010,10),DCB1-DCB)
plt.show()



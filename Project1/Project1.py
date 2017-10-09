# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 00:15:00 2017

@author: roco3
"""

import numpy as np
import math
# import matplotlib.pyplot as plt


# initial stock price
S0 = 90.74
# risk-free rate
r = 0.01
# dividend yield
div = 0.02
# sigma
sigma = 0.21
# down side threshold level
K = 72.592
# redemption price
red = [95.28, 99.81, 104.35]
# number of periods for a day
N = 750
# stock tree
VStock = np.zeros((1000,1000))
# price tree
VOption = np.zeros((1000,1000))
# barrie date
barr = 1/4

delta = 3/N

u = math.exp((r-div)*delta + sigma*math.sqrt(delta))
d = math.exp((r-div)*delta - sigma*math.sqrt(delta))

qu = (math.exp(r*delta)-d) / (u-d)
qd = 1 - qu

for i in range(N+1):
    VStock[N,i] = S0 * u**i * d**(N-i)
    if VStock[N,i] >= K:
        VOption[N,i] = 10.2125
    else:
        VOption[N,i] = 10 * VStock[N+1,i] / S0
for j in range(N-1, -1, -1):
    for i in range(j,-1,-1):
        VStock[j,i] = S0 * u**i * d**(j-i)
        VOption[j,i] = qu * VOption[j+1,i+1] + qd * VOption[j+1,i+1]
        for k in range(3):
            if j in k + np.array([1/4, 2/4,3/4,1]):
                if VStock[j,i] > K and VStock[j,i] <= red[k] :
                    VOption[j,i] = VOption[j,i] + 0.2125
                elif VStock[j,i] > red[k]:
                    VOption[j,i] = 10.2125

print(VOption[0,0])
        
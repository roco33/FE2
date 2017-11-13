# -*- coding: utf-8 -*-
"""
Created on Sun Nov 12 21:55:23 2017

@author: roco33
"""

import matplotlib.pyplot as plt
import numpy as np

s0 = 109.893


def VT(s):
    k = 1000 / s0
    D = 15.875
    if s < 0.7 * s0:
        V = k * s
    elif s >= 0.7 * s0 and s < s0:
        V = k * s + D
    else:
        V = 1000 + D
    return V


def VNT(s):
    D = 15.875
    if s < 0.7 * s0:
        V = 1000
    elif s >= 0.7 * s0:
        V = 1000 + D
    return V


vt = [VT(s) for s in range(130)]

vnt = [VNT(s) for s in range(130)]

plt.subplot(121)
plt.plot(np.arange(130), vt)
axes = plt.gca()
axes.set_ylim([0,1100])
axes.set_xlim([0,130])
plt.xlabel('Stock price')
plt.ylabel('Proceed to maturity')
plt.title('Triggered')

plt.subplot(122)
plt.plot(np.arange(130),vnt)
axes = plt.gca()
axes.set_ylim([0,1100])
axes.set_xlim([0,130])
plt.xlabel('Stock price')
# plt.ylabel('Proceed to maturity')
plt.title('Not triggered')
plt.show()

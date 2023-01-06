#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 27 00:09:05 2022

@author: kelvinleong
"""
import numpy as np
import matplotlib.pyplot as plt
from random import random
from math import exp, sqrt
#%%
def f(x):
    return x**(-1/2)/(1+exp(x))

a = 0
b = 1
N = 10000
################################### Q3a ######################################
MVM_I_arr = np.empty(100)
for k in range(100):
    fx_points = np.empty(N)
    for i in range(N):
        x_points = random()
        fx_points[i] = f(x_points)
    MVM_I_arr[k] = (b-a)/N * np.sum(fx_points)
    
MVM_I = np.mean(MVM_I_arr)
print(MVM_I)

################################### Q3b ######################################
#%%
def p(x):
    return 1/(2*sqrt(x))

def g(x, f, p):
    return f(x)/p(x)

IS_I_arr = np.empty(100)
for k in range(100):
    gx_points = np.empty(N)
    for i in range(N):
        x_points = random()**2
        gx_points[i] = g(x_points, f, p)
    IS_I_arr[k] = 1/N * np.sum(gx_points)
    
IS_I = np.mean(IS_I_arr)
print(IS_I)

################################## Q3c #######################################
#%%
bins = np.linspace(0.80, 0.8888, 10)

plt.figure(figsize=(8,8))
plt.hist(MVM_I_arr, bins, align='left')
#plt.xticks(bins[0:len(bins)-1])
plt.xlabel("Value",fontsize=16)
plt.ylabel("Distribution",fontsize=16)
plt.title("Mean value method MC integration",fontsize=16)

plt.figure(figsize=(8,8))
plt.hist(IS_I_arr, bins, align='left')
#plt.xticks(bins[0:len(bins)-1])
plt.xlabel("Value",fontsize=16)
plt.ylabel("Distribution",fontsize=16)
plt.title("Importance Sampling MC integration",fontsize=16)
plt.show()
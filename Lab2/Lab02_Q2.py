#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 22:56:22 2022

@Collaborator:  Kelvin Leong
                Patrick Sandoval
"""
####################### HEADER ###########################################

import numpy as np
import matplotlib.pyplot as plt
from time import time

pi = np.pi


# Define integrand as a function f(x)
def f(x):
    return 4/(1+x**2)

# Calculate area under the curve using trapezoidal rules (eq. 5.3)
def trapezoidal(a,b,N):
    # Width of slice
    h = (b - a) / N
    
    result = 0.5 * (f(a) + f(b))
    for k in range(1, N):
        result += f(a + k*h)
    result = h * result
    return result

# Calculate area under the curve using simpson's rules (eq. 5.9)
def simpsons(a,b,N):
    # Width of slice
    h = (b - a) / N
    
    odd_sum = 0
    even_sum = 0
    
    for k in range(1, N):
        if (k % 2)== 1:     # k is odd
            odd_sum += f(a + k*h)
        else:
            even_sum += f(a + k*h)
    result = h/3 * (f(a) + f(b) + 4 * odd_sum + 2 * even_sum)
    return result

#%%
################################ Q2b ####################################
N = 4
a = 0
b = 1

trapezoidal_area = trapezoidal(a,b,N)
simpsons_area = simpsons(a,b,N)

print(f"Exact value: {pi: .10f}")
print(f"Trapezoidal approximation: {trapezoidal_area: .10f}")
print(f"Simpsons approximation: {simpsons_area: .10f}")

#%%
################################ Q2c ####################################

n_simpsons = 0
a = 0
b = 1

for n in range(20):
    N = 2**n
    trapezoidal_area = trapezoidal(a,b,N)
    print(f"n={n} trapezoidal apprxoimation: {trapezoidal_area: .10f}")
print(f"Exact value: {pi: .10f}\n")



for n in range(20):
    N = 2**n
    simpsons_area = simpsons(a,b,N)
    print(f"n={n} simpsons apprxoimation: {simpsons_area: .10f}")
print(f"Exact value: {pi: .10f}\n")

# Manually finding
n_trapezoidal = 14
n_simpsons = 5
N_trapezoidal = 2**n_trapezoidal
N_simpsons = 2**n_simpsons

print(f"Number of slices for trapezoidal method to apprxomate integral with an error of O(10^-9): N={N_trapezoidal} slices (n={n_trapezoidal})")
print(f"Number of slices for simpsons method to apprxomate integral with an error of O(10^-9): N={N_simpsons} slices (n={n_simpsons})")

#%%
time_trapezoidal = []
time_simpsons = []

for i in range(500):
    start = time()
    _ = trapezoidal(a,b,N_trapezoidal)
    end = time()
    time_trapezoidal.append(end-start)
    
    start = time()
    _ = simpsons(a,b,N_simpsons)
    end = time()
    time_simpsons.append(end-start)
    
avg_time_trapezoidal = np.average(np.array(time_trapezoidal))
print(f"Average time to compute integral with trapezoidal method with an error of O(10^-9): {avg_time_trapezoidal:.5f} seconds")
    
avg_time_simpsons = np.average(np.array(time_simpsons))
print(f"Average time to compute integral with simpsons method with an error of O(10^-9): {avg_time_simpsons:.5f} seconds")


#%%
################################## Q2d #####################################
N1 = 16
N2 = 32
a = 0
b = 1

I1 = trapezoidal(a,b,N1)
I2 = trapezoidal(a,b,N2)

# Evaluating "practical estimation of errors" for N2 slice trapezoidal method (using eq.5.28)
error_N2_trapezoidal = (I2 - I1)/3
print(f"The practical error estimation trapezoidal method of N_2 = 32 slice is {error_N2_trapezoidal: .6f}")




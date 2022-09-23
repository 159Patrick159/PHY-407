#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 19 22:56:22 2022

@Collaborator:  Kelvin Leong
                Patrick Sandoval
"""
############################ HEADER ##########################################
# This python code contains the pseudo-code, codes and plots for calculating
# Trapezoidal and Simpson's approximation for the integral in the Q2 exercise
##############################################################################
# See pseudo-code in each of the sub-section

# Import the necessary libraries
import numpy as np
from time import time
import MyFunctions as myf

# From Q2a, the exact value of the integral is pi
pi = np.pi

#%%
################################ Q2b ####################################
# Pseudocode:
# Initialize conditions for integral: Lower bound=0, Upper bound=1, 4 slices for approximations
# Define the integrand as a function f(x)
# Calculate area under the curve using the Trapezoidal rule (eq.5.3)
# Calculate area under the curve using the Simpson's rule (eq.5.9)
# Print out the exact value, the Trapezoidal and the Simpson's approximation

N = 4
a = 0
b = 1

f = myf.q2_f
trapezoidal_area = myf.trapezoidal(f,a,b,N)
simpsons_area = myf.simpsons(f,a,b,N)

print(f"Exact value: {pi: .10f}")
print(f"Trapezoidal approximation: {trapezoidal_area: .10f}")
print(f"Simpsons approximation: {simpsons_area: .10f}")

#%%
################################ Q2c ####################################
# Pseudocode:
# (1) Find the number of slices to approximate integral with error O(10^-9)
# Run a for-loop with n-value between 1 and 20:
    # Let N = 2^n
    # Calculate area under the curve using the Trapezoidal rule (eq.5.3) with N slices
    # Print out the Trapezoidal area with N slices
# Run another for-loop with n-value between 1 and 20:
    # Let N = 2^n
    # Calculate area under the curve using the Simpson's rule (eq.5.9) with N slices
    # Print out the Simpson's area with N slices
# Print out the exact value of the integral, and manually find the n-value for the 
# Trapezoidal and Simpson's approximation in the printouts that matches the exact value
# until O(10^-9)
# Store these two values into n_trapezoidal and n_simpsons

# (2) Find average time it takes to compute the integral with an error of O(10−9)
# Run a for-loop with 500 iterations:
    # Begin the timing
    # Calculate trapezoidal area with N=2^(n_trapezoidal) slices
    # End the timing
    # Store the difference between end-time and start-time in time_trapezoidal array
    # Begin the timing
    # Calculate simpson's area with N=2^(n_simpsons) slices
    # End the timing
    # Store the difference between end-time and start-time in time_simpsons array
# Take the average of time_trapezoidal and time_simpsons arrays and print out results

for n in range(20):
    N = 2**n
    trapezoidal_area = myf.trapezoidal(f,a,b,N)
    print(f"n={n} trapezoidal apprxoimation: {trapezoidal_area: .10f}")
print(f"Exact value: {pi: .10f}\n")

for n in range(20):
    N = 2**n
    simpsons_area = myf.simpsons(f,a,b,N)
    print(f"n={n} simpsons apprxoimation: {simpsons_area: .10f}")
print(f"Exact value: {pi: .10f}\n")

# Manually finding the value
n_trapezoidal = 14
n_simpsons = 5
N_trapezoidal = 2**n_trapezoidal
N_simpsons = 2**n_simpsons

print(f"Number of slices for trapezoidal method to apprxomate integral with an error of O(10^-9): N={N_trapezoidal} slices (n={n_trapezoidal})")
print(f"Number of slices for simpsons method to apprxomate integral with an error of O(10^-9): N={N_simpsons} slices (n={n_simpsons})")


time_trapezoidal = []
time_simpsons = []

for i in range(500):
    start = time()
    _ = myf.trapezoidal(f,a,b,N_trapezoidal)
    end = time()
    time_trapezoidal.append(end-start)
    
    start = time()
    _ = myf.simpsons(f,a,b,N_simpsons)
    end = time()
    time_simpsons.append(end-start)
    
avg_time_trapezoidal = np.average(np.array(time_trapezoidal))
print(f"Average time to compute integral with trapezoidal method with an error of O(10^-9): {avg_time_trapezoidal*100:.5f} ms")
    
avg_time_simpsons = np.average(np.array(time_simpsons))
print(f"Average time to compute integral with simpsons method with an error of O(10^-9): {avg_time_simpsons*100:.5f} ms")


#%%
################################## Q2d #####################################
# Pseudocode:
# Define number of slices (N_1=16 slices, N_2=32 slices)
# Calculate the Trapezoidal approximation I_1 and I_2 using N_1 and N_2 slices respectively
# Evaluating "practical estimation of error" for N2 slice trapezoidal method (using eq.5.28)
# Print out the "practical estimation of error"

N1 = 16
N2 = 32

I1 = myf.trapezoidal(f,a,b,N1)
I2 = myf.trapezoidal(f,a,b,N2)

error_N2_trapezoidal = (I2 - I1)/3
print(f"The practical error estimation trapezoidal method of N_2 = 32 slice is {error_N2_trapezoidal: .6f}")




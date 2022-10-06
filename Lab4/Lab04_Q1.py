#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Contributors:
# Patrick Sandoval
# Kelvin Leong

################################ HEADER ######################################
# This file contains the main body code, pseudocode, and plots for solving the
# linear systems with different methods in Q1
##############################################################################
# Pseudo-code:
# Import the necessary libraries and functions from SolveLinear
# Implement the PartiatPivot() function according to textbook pg.222 
# (see SolveLinear.py for more detail)
# Q1a:
# Write matrix A and vector v from eq.6.2 in python
# Use GaussElim() and PartialPivot() with A and v to check if x gives correct output
# Q1b:
# Initialize an array of N from 5 to 250
# Initialize arrays of timing and arrays of error for each GaussElim, PartialPivot, LU method
# For each N in the N-array:
    # Randomly create NxN matrix A and N size vector v
    # Start timer
    # Run GaussElim/PartialPivot/LU to solve for x
    # End timer
    # Calculate the solution error by comparing v with vsol obtained by numpy.dot() with A and x
    # Append timings and error to arrays
# Plot the timings and errors for each N for each method

# Import libraries
import numpy as np
from numpy.random import rand 
from SolveLinear import GaussElim, PartialPivot
from numpy.linalg import solve
import matplotlib.pyplot as plt
from time import time

############################### Q1a ##########################################
# Please refer in the SolveLinear.py file for the implementation of the
# PartialPivot() function,

A = np.matrix([[2,1,4,1],[3,4,-1,-1],[1,-4,1,5],[2,-2,1,3]], float)
v = np.array([-4,3,9,7],float)
#A = np.matrix([[0,1,4,1],[3,4,-1,-1],[1,-4,1,5],[2,-2,1,3]], float)
#v = np.array([-4,3,9,7],float)
x = GaussElim(A,v)
print("Solution: ",x)
x = PartialPivot(A,v)
print("Solution: ",x)

#%%
############################### Q1b ##########################################
N_array = np.arange(5,251,1)

# Initialize arrays for plotting
time_GE = []
time_PP = []
time_LU = []
err_GE = []
err_PP = []
err_LU = []

for N in N_array:
    # Randomly generating matrix A and vector v
    A = rand(N,N)
    v = rand(N)

    # Timing how it takes to solve in Gaussian Elimination and corresponding error
    start = time()
    x = GaussElim(A,v)
    end = time()
    time_GE.append(end-start)
    vsol = np.dot(A, x)
    err = np.mean(abs(v-vsol))
    err_GE.append(err)
    
    # Timing how it takes to solve in Partial Pivoting and corresponding error
    start = time()
    x = PartialPivot(A,v)
    end = time()
    time_PP.append(end-start)
    vsol = np.dot(A, x)
    err = np.mean(abs(v-vsol))
    err_PP.append(err)
    
    # Timing how it takes to solve in LU decomposition and corresponding error
    start = time()
    x = solve(A,v)
    end = time()
    time_LU.append(end-start)
    vsol = np.dot(A, x)
    err = np.mean(abs(v-vsol))
    err_LU.append(err)
    
#%%
# Plotting timing as a function of N for each method
plt.figure(figsize=(8,8))
plt.plot(N_array, time_GE, 'ko', markersize=5, label='GaussElim')
plt.plot(N_array, time_PP, 'mo', markersize=5, label='PartialPivot')
plt.plot(N_array, time_LU, 'ro', markersize=5, label='LU')
plt.xlabel("N", fontsize=16)
plt.ylabel("Time [s]", fontsize=16)
plt.yscale('log')
plt.grid('on')
plt.legend()
plt.title("Time to solve N size linear system for each methods", fontsize=16)
plt.savefig('Q1_time.png')

# Plotting error as a function of N for each method
plt.figure(figsize=(8,8))
plt.plot(N_array, err_GE, 'ko', markersize=5, label='GaussElim')
plt.plot(N_array, err_PP, 'mo', markersize=5, label='PartialPivot')
plt.plot(N_array, err_LU, 'ro', markersize=5, label='LU')
plt.xlabel("N", fontsize=16)
plt.ylabel("Error", fontsize=16)
plt.yscale('log')
plt.grid('on')
plt.legend()
plt.title("Error in solving N size linear system for each methods", fontsize=16)
plt.savefig('Q1_error.png')

plt.show()
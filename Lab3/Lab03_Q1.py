#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code and pseudo-code for Q1 and all its subsections 
# as well as the plots and numerical values of our differentiation schemes.

######################################## Q1.a ############################################
# Pseudo-code:
# Import the necessary libraries and functions from MyFunctions.py
# Define array for h starting from 1e-16 to 1 as a power of 10
# Derive the analytic first derivative, and compute its value of x=0.5
# Calculate the numerical derivatives at x=0.5 with each value of h using forward 
# difference scheme from Eq.1 of lab manual (or Eq. 5.90 in textbook)
# Take the absolute difference between the analytic value and numerical derivatives as errors
# Plot the absolute error as a function of h on log-log plot

# Import functions from MyFunctions.py
from MyFunctions import erf, first_derivative_erf, fwd_derivative, central_derivative 
import numpy as np
import matplotlib.pyplot as plt

x_0 = 0.5
h = np.array([10**val for val in range(-16,1)])
#print(h)

analytic_val = first_derivative_erf(x_0)
fwd_numerical_derivatives = fwd_derivative(erf, x_0, h)
abs_fwd_numerical_error = abs(analytic_val - fwd_numerical_derivatives)

# Prints out analytic and numerical values
print(f"Analytic value of first derivative at x=0.5: f'(x)={analytic_val:.5f}")
print("Forward difference numerical derivative and error at x=0.5 for each h: ")
for i in range(len(h)):
    print(f"h = {h[i]}: f'(x){fwd_numerical_derivatives[i]:.5f}, error={abs_fwd_numerical_error[i]:.3e}")

####################################### Q1.b ##############################################
# Plotting the (absolute) numerical error as a function of h
plt.figure(figsize=(8,8))
plt.plot(h, abs_fwd_numerical_error, 'ro', markersize=5)
plt.xlabel('h', fontsize=16)
plt.ylabel('Absolute Numerical Error', fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.title("Numerical error of forward difference method at x=0.5", fontsize=16)
plt.grid('on')
plt.savefig('Q1b.png')

####################################### Q1.c ##############################################
#%%
# Pseudo-code:
# (Import the necessary functions from MyFunctions.py)
# (Define the same array of h as in Q1.a)
# (Derive the analytic first derivative, and compute its value of x=0.5)
# Calculate the numerical derivatives at x=0.5 with each value of h using central 
# difference scheme from Eq.2 of lab manual (or Eq. 5.98 in textbook)
# Take the absolute difference between the analytic value and numerical derivatives as errors
# Plot the absolute error of central difference scheme and error of forward difference scheme
# as a function of h on log-log plot

central_numerical_derivatives = central_derivative(erf, x_0, h)
abs_central_numerical_error = abs(analytic_val - central_numerical_derivatives)

# Plotting
plt.figure(figsize=(8,8))
plt.plot(h, abs_fwd_numerical_error, 'ro', markersize=5, label="Forward diff")
plt.plot(h, abs_central_numerical_error, 'bo', markersize=5, label="Central diff")
plt.xlabel('h', fontsize=16)
plt.ylabel('Absolute Numerical Error', fontsize=16)
plt.xscale('log')
plt.yscale('log')
plt.title("Numerical error of forward and central difference at x=0.5", fontsize=16)
plt.grid('on')
plt.legend(loc='lower left')
plt.savefig('Q1c.png')

plt.show()
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Contributors:
# Patrick Sandoval
# Kelvin Leong

################################ HEADER ######################################
# This file contains the main body code, pseudocode, and plots for solving the
# non-linear systems with different methods in Q3
##############################################################################
# See pseudo-code for each subsection
# Import libraries and functions
import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import Q3ab_f, Q3ab_f_derivative, relaxation_method, overrelaxation_method, Q3c_f, Q3c_f_derivative, binary_search, Newton_method
#%%
################################# Q3a ########################################
# Exercise 6.10 (a, b)
# Part a: Set c = 2, accuracy to 10^(-6)
# Gives a initial guess of x
# Write f(x) = 1 - e^{-cx} and its derivative as a function
# Implement the relaxation method:
    # While the error is greater than desired accuracy,
    # Record the value of x_next = f(x_guess) and x_guess,
    # Calculate the error using these two values and f'(x_next) according to Eq 6.83
    # Repeat until error less than desired accuracy
# Prints out the numerical solution x and its error
# Part b: Initialze an array of c from 0 to 3 in steps of 0.01
# Use relaxation method to store the numerical solution x at each c
# Plots of x as a function of c
    
# Part a
c = 2
accuracy = 1e-6
x_guess = 0.5

x, error, _ = relaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy)

print(f"Initial guess: {x_guess}")
print(f"Numerical solution for x: {x:.6f} converge to a solution accurate to 1e-6")
print(f"Error of the numerical solution: {error:.3e}")

# Part b
c_array = np.arange(0.01, 3, 0.01)
x_sol_array = []

for c in c_array:
    x_guess = 0.5
    x, error, _ = relaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy)

    x_sol_array.append(x)
    
# Plotting x as a function of c
plt.figure(figsize=(8,8))
plt.plot(c_array, x_sol_array, 'ro')
plt.xlabel("c", fontsize=16)
plt.ylabel("x", fontsize=16)
plt.grid('on')
plt.title('Numerical solution of x for 0 < c < 3',fontsize=16)
plt.savefig('Q3a.png')

#%%
################################ Q3b #########################################
# Exercise 6.11 (b,c,d)
# Pseudo-code: 
# Part b: Again set c = 2 and accuracy to 10^{-6}
# Use the implemented method to printout the numerical solution x,
# the error and this time also the counter for how many steps it takes to solve
# Part c: Implement of the over-relaxation method similar to relaxation method 
# but takes in omega as additional argument
# Modify x_next = (1+omega)*f(x_guess) - omega*x
# Modify error function as in part a
# Starts with omega = 0.5, manually change omega until counter is minimized
# Print out the numerical solution, error, and counter of steps
    
# Part b
c = 2
accuracy = 1e-6
x_guess = 0.5

x, error, counter = relaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy)

print(f"Initial guess: {x_guess}")
print(f"Numerical solution for x: {x:.6f}, converge to a solution accurate to 1e-6 after {counter} steps")
print(f"Error of the numerical solution: {error:.3e}")


# Part c
c = 2
accuracy = 1e-6
x_guess = 0.5
omega = 0.68

x, error, counter = overrelaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy, omega)

print("Over-relaxation iteration")
print(f"Initial guess: {x_guess}")
print(f"Numerical solution for x: {x:.6f}, converge to a solution accurate to 1e-6 after {counter} steps")
print(f"Error of the numerical solution: {error:.3e}")

#%%
############################### Q3c ##########################################
# Exercise 6.13
# Set accuracy to 10^{-6}
# Define function f(x) = 5e^{-x} + x - 5 (x such that f(x)=0 is a solution)
# and calculate its derivative function f'(x)
# Gives two initial guesses x1,x2 that bound the solution, check if 
# f(x1) and f(x2) have opposite sign
# Implement the binary search method (see MyFunctions for detail)
# Implement the newton's method (see MyFunctions for detail)
# Use binary search, newton's method and relaxation method with same initial guess to numerically solve
# for x, and print out their answer and number of steps to compare
# Use binary search's answer to calculate Wien's constant
# Use Wien's constant and Sun's peak wavelength to find temperature

accuracy = 1e-6

# Binary search method (procedure in pg.264)
x1 = 3.0
x2 = 5.0
print(f"f(x1) and f(x2) has opposite sign: {Q3c_f(x1)*Q3c_f(x2) < 0}")

binary_solution, binary_counter = binary_search(Q3c_f,x1,x2,accuracy, 0)
print(f"x solution from binaray search: {binary_solution:.7f}, in {binary_counter} steps")


### Q3d ##########
# Newton's method
newton_solution, newton_counter = Newton_method(Q3c_f, Q3c_f_derivative, x1, accuracy)
print(f"x solution from Newton's method: {newton_solution:.7f}, in {newton_counter} steps")

# Relaxation method
error = 1.0
counter = 0
while error > accuracy:
    x1, x2 = 5*(1-np.exp(-x1)), x1
    error = abs((x1 - x2)/(1- 1/5*np.exp(-x1)))
    #print(x1, error)
    counter += 1

relaxation_solution, relaxation_counter = x1, counter
print(f"x solution from relaxation method: {relaxation_solution:.7f}, in {relaxation_counter} steps")
###################


# Calculating Wien's constant
k_B = 1.380658e-23      # [m^2 kg s^-2 K^-1], Boltzmann Constant
h = 6.6260755e-34       # [m^2 kg s], Planck's constant
c = 2.99792458e8        # [m/s]
x = binary_solution

b = (h*c)/(k_B*x)       # [m K]
print(f"Numerically calculated Wien's constant is b = {b: .7f} m*K")

Sun_pk_wvlength = 502e-9 #[m]
T_sun = b/Sun_pk_wvlength
print(f"Surface temperature of the Sun is T = {T_sun: .3f} K")



#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Contributors:
# Patrick Sandoval
# Kelvin Leong

################################ HEADER ######################################
# This file contains the main body code, pseudocode, and plots for solving the
# non-linear systems with different methods in Q3
##############################################################################

import numpy as np
import matplotlib.pyplot as plt

################################# Q3a ########################################
# Exercise 6.10 (a, b)
# Part a
c = 2
accuracy = 1e-6
x_guess = 0.5
error = 1.0

print(f"Initial guess: {x_guess}")

def Q3ab_f(x,c):
    return 1 - np.exp(-c*x)

def Q3ab_f_derivative(x,c):
    return c*np.exp(-c*x)

def relaxation_method(f, fprime, x, c, accuracy):
    x1 = x
    error = 1.0
    counter = 0
    # Loop until error is smaller than target accuracy
    while error > accuracy:
        x1, x2 = f(x1,c), x1
        error = abs((x1 - x2)/(1- 1/fprime(x1,c)))
        #print(x1, error)
        counter += 1
    return x1, error, counter

x, error, _ = relaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy)
#while error > accuracy:
#    x1, x2 = 1 - np.exp(-c*x1), x1
#    error = abs((x1 - x2)/(1- 1/(c*np.exp(-c*x2))))
    #print(x1, error)

print(f"Numerical solution for x: {x:.6f} converge to a solution accurate to 1e-6")
print(f"Error of the numerical solution: {error:.3e}")

#%%
# Part b
c_array = np.arange(0.01, 3, 0.01)
x_sol_array = []

for c in c_array:
    x_guess = 0.5
    x, error, _ = relaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy)
    #error = 1.0
    # Loop until error is smaller than target accuracy
    #while error > accuracy:
    #    x1, x2 = 1 - np.exp(-c*x1), x1
    #    error = abs((x1 - x2)/(1- 1/(c*np.exp(-c*x2))))
    #    #print(x1, error)
    x_sol_array.append(x)
    

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
# Part b
c = 2
accuracy = 1e-6
x_guess = 0.5
#error = 1.0
#steps = 0

print(f"Initial guess: {x_guess}")
x, error, counter = relaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy)
#while error > accuracy:
#    x1, x2 = 1 - np.exp(-c*x1), x1
#    error = abs((x1 - x2)/(1- 1/(c*np.exp(-c*x2))))
    #print(x1, error)
#    steps += 1
print(f"Numerical solution for x: {x:.6f}, converge to a solution accurate to 1e-6 after {counter} steps")
print(f"Error of the numerical solution: {error:.3e}")

#%%
# Part c
c = 2
accuracy = 1e-6
x_guess = 0.5
#error = 1.0
omega = 0.68
#steps = 0

def overrelaxation_method(f, fprime, x, c, accuracy, omega):
    x1 = x
    error = 1.0
    counter = 0
    # Loop until error is smaller than target accuracy
    while error > accuracy:
        x1, x2 = (1+omega)*f(x1,c) - omega*x1, x1
        error = abs((x1 - x2)/(1- 1/((1+omega)*fprime(x1,c) - omega)))
        #print(x1, error)
        counter += 1
    return x1, error, counter

x, error, counter = overrelaxation_method(Q3ab_f, Q3ab_f_derivative, x_guess, c, accuracy, omega)

print("Over-relaxation iteration")
print(f"Initial guess: {x_guess}")
#while error > accuracy:
#    x1, x2 = (1+omega)*(1 - np.exp(-c*x1))-omega*x1, x1
#    error = abs((x1 - x2)/(1- 1/((1+omega) * c*np.exp(-c*x2) - omega)))
#    steps += 1
print(f"Numerical solution for x: {x:.6f}, converge to a solution accurate to 1e-6 after {counter} steps")
print(f"Error of the numerical solution: {error:.3e}")

#%%
############################### Q3c ##########################################
# Exercise 6.13
# x is the solution to f(x) = 5e^{-x} + x - 5 = 0 
# 
accuracy = 1e-6

# Binary search method (procedure in pg.264)
x1 = 3.0
x2 = 5.0

def Q3c_f(x):
    return 5*np.exp(-x) + x - 5

def Q3c_f_derivative(x):
    return -5*np.exp(-x) + 1

# Check if f(x1) and f(x2) have opposite signs
print(f"f(x1) and f(x2) has opposite sign: {Q3c_f(x1)*Q3c_f(x2) < 0}")

# Use recursion to implement binary search
def binary_search(f, x1, x2, accuracy):
    global binary_counter
    binary_counter += 1
    
    f_x1 = f(x1)
    #f_x2 = f(x2)
    x_half = 0.5 * (x1 + x2)
    f_x_half = f(x_half)
    
    if (f_x1 * f_x_half > 0):
        x1 = x_half
    else:
        x2 = x_half
        
    if(abs(x1 - x2) < accuracy):
        return 0.5*(x1 + x2)
    else:
        return binary_search(f, x1, x2, accuracy)
    
binary_counter = 0
binary_solution = binary_search(Q3c_f,x1,x2,accuracy)
print(f"x solution from binaray search: {binary_solution:.7f}, in {binary_counter} steps")


# Newton's method
def Newton_method(f, fprime, x, accuracy):
    counter = 0
    delta = 1.0
    while abs(delta) > accuracy:
        delta = f(x)/fprime(x)
        x -= delta
        counter +=1
    return x, counter

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



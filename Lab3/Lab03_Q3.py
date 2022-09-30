#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Contributors:
# Patrick Sandoval
# Kelvin Leong

################################# HEADER #####################################
# This python script contains holds the code and pseudo-code for Q3 and all its
# subsections, as well as the plots of wavefunctions of the quautum SHM
##############################################################################
# Pseudo-code:
# Import the necessary libraries and functions from MyFunctions
# Define the Hermite Polynomial as a function of (n,x) from Eq.9
# Define the wavefunction using the Hermite Polynomial from Eq.8
# (Q3.a) Define the domain of the wavefunction -4 < x < 4 (with 100 points)
# Get the wavefunction for n = {0,1,2,3} on the domain, and plot it on the same graph
# (Q3.b) Define the domain of the wavefunction -10 < x < 10 (with 1000 points)
# Get the wavefunction for n = 30 on this domain, and plot it on a separate graph
# (Q3.c)
# Define the first derivative of wavefunction from Eq.11 using the the H_n 
# and H_{n-1} polynomial
# To use evaluate the double improper integrals in Eq.12-13 for each n:
    # We first use Gaussian quadrature using gaussxwab with N=100 points and set 
    # the bounds a=-pi/2 and b=pi/2 (according to the evaluation technique from 
    # Eq 5.74-5.75) to obtain the new domain and weights
    # Then we define the after change of variables as functions of n and the new 
    # domain according to the evaluation technique from Eq. 5.75 in the textbooks
    # For n in {0,1,...,15}:
        # Evaluate the integrals using the integrand functions on the domain 
        # and with weights similar to the example in pg.180 in the textbook
# Print out the results of the integrals as <x^2> and <p^2> for each n
# Take the square root of results as their uncertainties
# Calculate the total energy of the system for each nth energy level with Eq. 14

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import gaussxwab, Hermite_Polynomial, Wave_function, \
    integrand_expected_position_square, integrand_expected_momentum_square


#%%
################################### Q3.a ##########################################
# Define domain -4 < x < 4
x_array = np.linspace(-4, 4, 100)

# Hermite Polynomials for n={0,1,2,3}
H0 = Hermite_Polynomial(0, x_array)
H1 = Hermite_Polynomial(1, x_array)
H2 = Hermite_Polynomial(2, x_array)
H3 = Hermite_Polynomial(3, x_array)
#H4 = Hermite_Polynomial(4, x_array)

# Wavefunctions for n={0,1,2,3}
Psi_0 = Wave_function(0, x_array)
Psi_1 = Wave_function(1, x_array)
Psi_2 = Wave_function(2, x_array)
Psi_3 = Wave_function(3, x_array)
#Psi_4 = Wave_function(4, x_array)

# Plotting the Hermite functions (not needed for submission)
plt.figure(figsize=(8,8))
plt.plot(x_array, H0, 'k--', label="H0")
plt.plot(x_array, H1, 'g--', label="H1")
plt.plot(x_array, H2, 'r--', label="H2")
plt.plot(x_array, H3, 'b--', label="H3")
#plt.plot(x_array, H4, 'm--', label="H4")
plt.xlabel("x", fontsize=16)
plt.ylabel("H(n,x)", fontsize=16)
plt.ylim([-50,50])
plt.title("Hermite Polynomial for n={0,1,2,3}",fontsize=16)
plt.grid('on')
plt.legend()

# Plotting the wavefunctions
plt.figure(figsize=(8,8))
plt.plot(x_array, Psi_0, 'k--', label=r"$\Psi_0(x)$")
plt.plot(x_array, Psi_1, 'g--', label=r"$\Psi_1(x)$")
plt.plot(x_array, Psi_2, 'r--', label=r"$\Psi_2(x)$")
plt.plot(x_array, Psi_3, 'b--', label=r"$\Psi_3(x)$")
#plt.plot(x_array, Psi_4, 'm--', label=r"$Psi_4(x)$")
plt.xlabel("x", fontsize=16)
plt.ylabel(r"$\Psi_n(x)$",fontsize=16)
plt.title(r"Wavefunction $\Psi_n(x)$ for n={0,1,2,3}",fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig('Q3a.png')

#%%
################################### Q3.b ##########################################
# Define domain -4 < x < 4
x_array = np.linspace(-10, 10, 1000)

# Wavefunction at n=30
Psi_30 = Wave_function(30, x_array)

# Plotting the wavefunctions at n=30
plt.figure(figsize=(8,8))
plt.plot(x_array, Psi_30, 'k-', label=r"$\Psi_{30}(x)$")
plt.xlabel("x", fontsize=16)
plt.ylabel(r"$\Psi_n(x)$",fontsize=16)
plt.title(r"Wavefunction $\Psi_n(x)$ for n=30",fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig('Q3b.png')

#%%
################################### Q3.c ##########################################
# Set n={0,...,15}
n_array = np.arange(0, 16, 1)

# Gaussian quadrature with N=100 (only need to calculate once)
N = 100
a = -np.pi/2
b = np.pi/2
x_array, w = gaussxwab(N, a, b)

# Arrays to store relavent results for Q3c
expected_x2_array = []
expected_p2_array = []
uncertainty_x2_array = []
uncertainty_p2_array = []
E_array = []

for n in n_array:
    expected_x2 = 0.0
    expected_p2 = 0.0
    # For-loop to evaluate integral at n
    for k in range(N):
        expected_x2 += w[k] * integrand_expected_position_square(n, x_array[k])
        expected_p2 += w[k] * integrand_expected_momentum_square(n, x_array[k])
    uncertainty_x2 = np.sqrt(expected_x2)
    uncertainty_p2 = np.sqrt(expected_p2)
    E_n = 0.5 * (expected_x2 + expected_p2)
    
    # Append the results to array
    expected_x2_array.append(expected_x2)
    expected_p2_array.append(expected_p2)
    uncertainty_x2_array.append(uncertainty_x2)
    uncertainty_p2_array.append(uncertainty_p2)
    E_array.append(E_n)
    
    print(f"n={n}: <x^2>={expected_x2:.3f}, sqrt(<x^2>)={uncertainty_x2:.3f}; <p^2>={expected_p2:.3f}, sqrt(<p^2>)={uncertainty_p2:.3f}; E={E_n:.3f}")


plt.show()



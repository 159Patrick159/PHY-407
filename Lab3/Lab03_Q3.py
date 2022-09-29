#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Contributors:
# Patrick Sandoval
# Kelvin Leong

################################# HEADER #####################################
# This python script contains holds the code and pseudo-code for Q3 and all its
# subsections
##############################################################################
import numpy as np
import matplotlib.pyplot as plt
from math import factorial
import MyFunctions as myf


def Hermite_Polynomial(n,x):
    """
    Parameters: n: nth energy level, integer greater equal than 0
                x: Position
    Return: nth Hermite Polynomial
    """
    H_0 = np.ones(len(x))
    H_1 = 2*x
    if n == 0: return H_0
    if n == 1: return H_1
    
    H_array = [H_0, H_1]
    
    for i in range(1, n):
        H_next = 2*x*H_array[i] - 2*i*H_array[i-1]
        H_array.append(H_next)
        
    return H_array[-1]

def Wave_function(n, x):
    H_n = Hermite_Polynomial(n, x)
    return np.exp(-x**2/2)/np.sqrt(2**n * factorial(n) * np.sqrt(np.pi)) * H_n

################################### Q3.a ##########################################
x_array = np.linspace(-4, 4, 100)
H0 = Hermite_Polynomial(0, x_array)
H1 = Hermite_Polynomial(1, x_array)
H2 = Hermite_Polynomial(2, x_array)
H3 = Hermite_Polynomial(3, x_array)
H4 = Hermite_Polynomial(4, x_array)

Psi_0 = Wave_function(0, x_array)
Psi_1 = Wave_function(1, x_array)
Psi_2 = Wave_function(2, x_array)
Psi_3 = Wave_function(3, x_array)
Psi_4 = Wave_function(4, x_array)

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

x_array = np.linspace(-10, 10, 1000)
Psi_30 = Wave_function(30, x_array)

plt.figure(figsize=(8,8))
plt.plot(x_array, Psi_30, 'k-', label=r"$\Psi_{30}(x)$")
plt.xlabel("x", fontsize=16)
plt.ylabel(r"$\Psi_n(x)$",fontsize=16)
plt.title(r"Wavefunction $\Psi_n(x)$ for n=30",fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig('Q3b.png')

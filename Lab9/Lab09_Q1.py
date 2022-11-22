#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 00:10:50 2022

@author: kelvinleong
"""

import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import calc_normalization_constant, V, psi

# Constants
h_bar = 1.055e-34   #[m^2 kg/s]
L = 1e-8            #[m]
m = 9.109e-31       #[kg]
sigma = L/25.0
kappa = 500/L
x0 = L/5

# Spatial quantities
P = 1024
a = L/P
p_arr = np.arange(1, P, 1)
x_p_arr = p_arr*a - L/2

# Time quantities
N = 3000
tau = 1e-18
T = N*tau

# Define potential V(x)
potential_arr = V(x_p_arr, -L/2, L/2)

# Define the discretized Hamiltonian matrix
A = -h_bar**2/(2*m*a**2)
B = potential_arr - 2*A
Hamiltonian = np.zeros([len(p_arr), len(p_arr)])
for i in range(len(p_arr)):
    if i == 0:
        Hamiltonian[i][i] = B[i]
        Hamiltonian[i][i+1] = A
    elif i == len(p_arr)-1:
        Hamiltonian[i][i-1] = A
        Hamiltonian[i][i] = B[i]
    else:
        Hamiltonian[i][i-1] = A
        Hamiltonian[i][i] = B[i]
        Hamiltonian[i][i+1] = A

# Define the initial (t=0) wave fucntion from psi_1 to psi_{P-1}
# and normalize it (using relation in eq3)
psi_0 = psi(x_p_arr, x0, sigma, kappa)
#nor_const = calc_normalization_constant(psi, P, -L/2, L/2, x0, sigma, kappa)
#psi_0_nor = psi_0 / np.sqrt(nor_const)

# Construct L and R matrices used for Crank-Nicolson time-stepper method
L_matrix = np.eye(P-1) + 1j * tau/(2*h_bar) * Hamiltonian.copy()
L_matrix_inv = np.linalg.inv(L_matrix)
R_matrix = np.eye(P-1) - 1j * tau/(2*h_bar) * Hamiltonian.copy()

# Define the receiver of solution, and give initial condition
psi_all = np.zeros((N, P-1), dtype=complex)      #[Time][Spatial]
psi_all[0] = psi_0

# Loop through time from n=0 to n=N-1 to solve for n+1 at each step
for n in range(N-1):
    psi_all[n+1, : ] = np.matmul(L_matrix_inv, np.matmul(R_matrix, psi_all[n, :]))

# End of solving the time-dependent Schrodinger equation
#%%
################################### Q1b ######################################
x_domain = x_p_arr
t_domain = np.linspace(0, T, N)

# Define indices of times for plotting
n0, n1, n2, n3, n4 = 0, int(N/4), int(N/2), int(3*N/4), N


plt.figure(figsize=(10,10))
plt.plot(x_domain, np.abs(psi_0))
plt.xlabel("Spatial Domain", fontsize=16)
plt.ylabel("Probability amplitude", fontsize=16)
plt.grid('on')

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 00:10:50 2022

@author: kelvinleong
"""

import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import calc_normalization_constant, V, psi

# Constants
h_bar = 1.055e-34   #[m^2 kg/s]
L = 1e-8            #[m]
m = 9.109e-31       #[kg]
sigma = L/25.0
kappa = 500/L
x0 = L/5

# Spatial quantities
P = 1024
a = L/P
p_arr = np.arange(1, P, 1)
x_p_arr = p_arr*a - L/2

# Time quantities
N = 3000
tau = 1e-18
T = N*tau

# Define potential V(x)
potential_arr = V(x_p_arr, -L/2, L/2)

# Define the discretized Hamiltonian matrix
A = -h_bar**2/(2*m*a**2)
B = potential_arr - 2*A
Hamiltonian = np.zeros([len(p_arr), len(p_arr)])
for i in range(len(p_arr)):
    if i == 0:
        Hamiltonian[i][i] = B[i]
        Hamiltonian[i][i+1] = A
    elif i == len(p_arr)-1:
        Hamiltonian[i][i-1] = A
        Hamiltonian[i][i] = B[i]
    else:
        Hamiltonian[i][i-1] = A
        Hamiltonian[i][i] = B[i]
        Hamiltonian[i][i+1] = A

# Define the initial (t=0) wave fucntion from psi_1 to psi_{P-1}
# and normalize it (using relation in eq3)
psi_0 = psi(x_p_arr, x0, sigma, kappa)
#nor_const = calc_normalization_constant(psi, P, -L/2, L/2, x0, sigma, kappa)
#psi_0_nor = psi_0 / np.sqrt(nor_const)

# Construct L and R matrices used for Crank-Nicolson time-stepper method
L_matrix = np.eye(P-1) + 1j * tau/(2*h_bar) * Hamiltonian.copy()
L_matrix_inv = np.linalg.inv(L_matrix)
R_matrix = np.eye(P-1) - 1j * tau/(2*h_bar) * Hamiltonian.copy()

# Define the receiver of solution, and give initial condition
psi_all = np.zeros((N, P-1), dtype=complex)      #[Time][Spatial]
psi_all[0] = psi_0

# Loop through time from n=0 to n=N-1 to solve for n+1 at each step
for n in range(N-1):
    psi_all[n+1, : ] = np.matmul(L_matrix_inv, np.matmul(R_matrix, psi_all[n, :]))

# End of solving the time-dependent Schrodinger equation
#%%
################################### Q1b ######################################
x_domain = x_p_arr
t_domain = np.linspace(0, T, N)

# Define indices of times for plotting
n0, n1, n2, n3, n4 = 0, int(N/4), int(N/2), int(3*N/4), N


plt.figure(figsize=(10,10))
plt.plot(x_domain, np.abs(psi_0))
plt.xlabel("Spatial Domain", fontsize=16)
plt.ylabel("Probability amplitude", fontsize=16)
plt.grid('on')


 
 
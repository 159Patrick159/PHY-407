#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Contributors 
# Patrick Sandoval
# Kelvin Leong

################################# HEADER ##########################################
# This python script holds the code, pseudo-code and plots for Q1 of Lab08 for    #  
# solving the potential inside a box with capacitor plates that satisfies the     #
# Laplace equations and boundary conditions                                       #
###################################################################################
# Pseudocode:
# Define a 101 x 101 grid points matrix, target accuracy
# Define the capacitor plate at x=2cm and x=8cm as a boundary condition
# Initialize delta=1.0
# while delta > target accuracy:
    # Loop through each grid point:
        # if grid point is at specified boundary
            # Does not change for next step
        # else:
            # Calculate average from adjacent value for next step (and with omega
            # if using over-relaxation)
        # Update delta as the maximum of difference in current step and next step
# Plot the contour of the potential solution and print how many iterations it takes to solve

import numpy as np
import matplotlib.pyplot as plt

#%%
################################# Q1a ########################################
# Constants
M = 100
V = 1.0     # Voltage un the capacitor 
target = 1e-6   # Target accuracy 

#phi, iteration = SolveCapacitor(V, target, 0, verbose=True)
phi = np.zeros([M+1, M+1], dtype=float)
phi[20:80, 20] = V
phi[20:80, 80] = -V

delta = 1.0
iteration = 0
verbose = True
err_array = np.zeros([M+1, M+1])

while delta > target:
    
    for i in range(M+1):
        for j in range(M+1):
            # Boundary Conditions
            if not (i == 0 or i == M or j == 0 or j == M or ((j == 20 or j == 80) and i in np.arange(20, 80, 1))):
                # Gauss-Seidel, over-relaxation method
                adjacent = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1]+ phi[i][j-1])
                
                old = phi[i][j]
                phi[i][j] = adjacent/4
                
                err_array[i][j] = abs(old - phi[i][j])
    
    # Calculate difference between old and new values
    delta = np.max(err_array)
    
    iteration += 1
    if verbose: print(f"Iteration {iteration}: delta = {delta}")
 


levels = np.arange(-1.0, 1.01, 0.02)
plt.figure(figsize=(12,10))
plt.contourf(np.arange(0, M+1, 1), np.arange(0, M+1, 1), phi, levels, cmap='coolwarm')
cbar = plt.colorbar()
cbar.set_label('Potential $V$', fontsize=20)
plt.xlabel('X [mm]', fontsize=20)
plt.ylabel('Y [mm]', fontsize=20)
plt.title('Contour map of potential $V$ of a capacitor', fontsize=20)
plt.axis('equal')
plt.tight_layout()
plt.savefig('Q1a.png')

################################# Q1c ########################################
Ey, Ex = np.gradient(-phi, np.arange(0, M+1, 1), np.arange(0, M+1, 1))
fig = plt.figure(figsize=(12,10))
strm = plt.streamplot(np.arange(0, M+1, 1), np.arange(0, M+1, 1), Ex, Ey, color=phi, linewidth=5, cmap='coolwarm')
cbar = fig.colorbar(strm.lines)
cbar.set_label('Potential $V$', fontsize=20)
plt.title('Electric field lines of a capacitor', fontsize=20)
plt.xlabel('X [mm]', fontsize=20)
plt.ylabel('Y [mm]', fontsize=20)
plt.axis('equal')
plt.tight_layout()
plt.savefig('Q1c.png')
#%%
################################## Q1b #######################################

#V = 1.0     # Voltage un the capacitor 
#target = 1e-6   # Target accuracy 

#phi, iteration = SolveCapacitor(V, target, 0.1, verbose=True)


#phi, iteration = SolveCapacitor(V, target, 0.5, verbose=True)
#print(f"Gauss-Seidel with over-relaxation method (omega = 0.5) took {iteration} iterations to solve")

#Create arrays to hold potential values
phi = np.zeros([M+1, M+1], dtype=float)
phi[20:80, 20] = V
phi[20:80, 80] = -V

delta = 1.0
iteration_1 = 0
verbose = True
err_array = np.zeros([M+1, M+1])

omega_1 = 0.1
while delta > target:
    
    for i in range(M+1):
        for j in range(M+1):
            # Boundary Conditions
            if not (i == 0 or i == M or j == 0 or j == M or ((j == 20 or j == 80) and i in np.arange(20, 80, 1))):
                # Gauss-Seidel, over-relaxation method
                adjacent = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1]+ phi[i][j-1])
                
                old = phi[i][j]
                phi[i][j] = adjacent * (1+omega_1)/4 - omega_1*old
                
                err_array[i][j] = abs(old - phi[i][j])
    
    # Calculate difference between old and new values
    delta = np.max(err_array)
    
    iteration_1 += 1
    if verbose: print(f"Iteration {iteration_1}: delta = {delta}")

#%%
levels = np.arange(-1.0, 1.01, 0.02)
plt.figure(figsize=(12,10))
plt.contourf(np.arange(0, M+1, 1), np.arange(0, M+1, 1), phi, levels, cmap='coolwarm')
cbar = plt.colorbar()
cbar.set_label('Potential $V$', fontsize=20)
plt.xlabel('X [mm]', fontsize=20)
plt.ylabel('Y [mm]', fontsize=20)
plt.title('Contour map of potential $V$ of a capacitor with omega=0.1', fontsize=20)
plt.axis('equal')
plt.tight_layout()
plt.savefig('Q1b-1.png')
#%%
#Create arrays to hold potential values
phi = np.zeros([M+1, M+1], dtype=float)
phi[20:80, 20] = V
phi[20:80, 80] = -V

delta = 1.0
iteration_5 = 0
verbose = True
err_array = np.zeros([M+1, M+1])

omega_5 = 0.5
while delta > target:
    
    for i in range(M+1):
        for j in range(M+1):
            # Boundary Conditions
            if not (i == 0 or i == M or j == 0 or j == M or ((j == 20 or j == 80) and i in np.arange(20, 80, 1))):
                # Gauss-Seidel, over-relaxation method
                adjacent = (phi[i+1][j] + phi[i-1][j] + phi[i][j+1]+ phi[i][j-1])
                
                old = phi[i][j]
                phi[i][j] = adjacent * (1+omega_5)/4 - omega_5*old
                
                err_array[i][j] = abs(old - phi[i][j])
    
    # Calculate difference between old and new values
    delta = np.max(err_array)
    
    iteration_5 += 1
    if verbose: print(f"Iteration {iteration_5}: delta = {delta}")


levels = np.arange(-1.0, 1.01, 0.02)
plt.figure(figsize=(12,10))
plt.contourf(np.arange(0, M+1, 1), np.arange(0, M+1, 1), phi, levels, cmap='coolwarm')
cbar = plt.colorbar()
cbar.set_label('Potential $V$', fontsize=20)
plt.xlabel('X [mm]', fontsize=20)
plt.ylabel('Y [mm]', fontsize=20)
plt.title('Contour map of potential $V$ of a capacitor with omega=0.5', fontsize=20)
plt.axis('equal')
plt.tight_layout()
plt.savefig('Q1b-2.png')

print(f"Gauss-Seidel with relaxation method took {iteration} iterations to solve")        
print(f"Gauss-Seidel with over-relaxation method (omega = {omega_1}) took {iteration_1} iterations to solve")
print(f"Gauss-Seidel with over-relaxation method (omega = {omega_5}) took {iteration_5} iterations to solve")
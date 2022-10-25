#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 10:21:28 2022

@author: kelvinleong
"""

import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import Verlet_multibody, KE, Lennard_Jones

#%%
# Number of particles simulating
N = 16

# Set initial conditions for the N particles
Lx = 4.0
Ly = 4.0

dx = Lx/np.sqrt(N)
dy = Ly/np.sqrt(N)

x_grid = np.arange(dx/2, Lx, dx)
y_grid = np.arange(dy/2, Ly, dy)

xx_grid, yy_grid = np.meshgrid(x_grid, y_grid)

x_initial = xx_grid.flatten()
y_initial = yy_grid.flatten()


# Assign the initial velocities of all N particles to zero
vx_initial = np.zeros(N)
vy_initial = np.zeros(N)

# Define time step 
dt = 0.01
T = 1000*dt
t = np.arange(0,T,dt)

r_initial = np.array([x_initial, y_initial])
v_initial = np.array([vx_initial, vy_initial])

x_all_t, y_all_t, vx_all_t, vy_all_t, E_t = Verlet_multibody(t,dt,r_initial,v_initial,N)

#%%
fig, (a0,a1,a2) = plt.subplots(figsize=(12,16),ncols=1,nrows=3,\
                               gridspec_kw={'height_ratios': [3, 1, 1]})

for i in range(len(x_all_t[0])):
    a0.plot(x_all_t.T[i],y_all_t.T[i],'.')
    a0.set_title("N=16 particles Simulation")
    a0.set_xlabel("X-Position")
    a0.set_ylabel("Y-Position")
    a0.grid('on')
    #a0.legend()
    
    a1.plot(t,x_all_t.T[i])
    a1.set_ylabel("X-Component")
    a1.set_xlabel("Time [s]")
    a1.grid('on')
    
    a2.plot(t,y_all_t.T[i])
    a2.set_ylabel("Y-Component")
    a2.set_xlabel("Time [s]")
    a2.grid('on')
    
plt.tight_layout()
plt.savefig("Q2a.png")

avg = np.median(E_t)
plt.figure(figsize=(8,8))
plt.plot(t,E_t, 'r.')
plt.plot(t, np.ones(len(t))*avg, 'k--', label='Mean Energy')
plt.plot(t, np.ones(len(t))*avg*1.01, 'm--', label='+/- 1% of mean')
plt.plot(t, np.ones(len(t))*avg*0.99, 'm--')
plt.xlabel('Time [s]',fontsize=16)
plt.ylabel('Energy',fontsize=16)
plt.title('N=16 particles system energy', fontsize=16)
plt.grid('on')
plt.legend()
plt.tight_layout()
plt.savefig("Q2b.png")
plt.show() 


    
                



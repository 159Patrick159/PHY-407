#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 11:04:32 2022

@Collaborator:  Kelvin Leong
                Patrick Sandoval
"""

#%%
############################### Q2a ##########################################
# Pseudocode:
    # Import the necessary libraries
    # Define all needed physical constants given in the physics background in lab
    # Define a time step dt, end_time
    # Create a time array t based on the defined dt and end_time
    # Initialize arrays of zeros x, y, vx, and vy with length the same as t for Earth and Jupiter respectively
    # Set the first element of the arrays to the initial conditions of Earth and Jupiter given in lab
    # (Use Euler-Cromer)
    # Run a for-loop from 0 to len(t)-1 iterations:
        # Compute current distance of Jupiter from Sun using x_J[i] and y_J[i]
        # Calculate current accelerations ax and ay based on eq(6) using current position x_J[i] and y_J[i]
        # Calculate velocities of next increment vx_J[i+1] and vy_J[i+1] using ax, ay, vx_J[i], vy_J[i], dt
        # Calculate positions of next increment x_J[i+1] and y_J[i+1] using vx_J[i], vy_J[i], x_J[i], y_J[i], dt
    # Run a for-loop from 0 to len(t)-1 iterations:
        # Compute current distance of Earth from Sun using x_E[i] and y_E[i]
        # Compute current distance of Earth from Jupiter using x_E[i], y_E[i], x_J[i], y_J[i]
        # Calculate current accelearations ax and ay due to Sun and Jupiter (note vector directions) 
        # using x_E[i], y_E[i], x_J[i], y_J[i]
        # Calculate velocities of next increment vx_E[i+1] and vy_E[i+1] using ax, ay, vx_E[i], vy_E[i], dt
        # Calculate positions of next increment x_E[i+1] and y_E[i+1] using vx_E[i], vy_E[i], x_E[i], y_E[i], dt
    # Plot the X v.s. Y position plot of Jupiter and Earth for 10 years on the graph

# Import libraries
import numpy as np
import matplotlib.pyplot as plt
import MyFunctions as myf

# Constants
G = 39.5        # [AU^3 * M_sol^-1 yr^-2] Graviational Constants 
M_sol = 1.0     # Solar Mass
M_J = 10e-3     # Jupiter Mass
alpha = 0.01    # [AU^2] Precession

# Function to calculate the distance between Earth and Jupiter
#def earth_jupiter_distance(x_E, y_E, x_J, y_J):
#    dx = np.abs(x_E - x_J)
#    dy = np.abs(y_E - y_J)
#    return np.sqrt(dx**2 + dy**2)


# Build time array
dt = 0.0001             # [year]
end_time = 10           # [year]
t = np.arange(0,end_time,dt)

# Initialize zero arrays for Earth and Jupiter
x_E = np.zeros(len(t))
y_E = np.zeros(len(t))
vx_E = np.zeros(len(t))
vy_E = np.zeros(len(t))

x_J = np.zeros(len(t))
y_J = np.zeros(len(t))
vx_J = np.zeros(len(t))
vy_J = np.zeros(len(t))

# Initial Condition for Earth and Jupiter
x_E[0] = 1.0    # [AU]
y_E[0] = 0.0    # [AU]
vx_E[0] = 0.0   # [AU/yr]
vy_E[0] = 6.18  # [AU/yr]
x_J[0] = 5.2    # [AU]
y_J[0] = 0.0    # [AU]
vx_J[0] = 0.0   # [AU/yr]
vy_J[0] = 2.63  # [AU/yr]

# Using Euler-Cromer method to simulate Jupiter's orbit due to Sun 
# (assuming gravitational force due to Earth is small) first
for i in range(0, len(t)-1):
    r_J = np.sqrt(x_J[i]**2 + y_J[i]**2)
    
    # Calculate acceleration
    ax = -G * M_sol * x_J[i]/r_J**3
    ay = -G * M_sol * y_J[i]/r_J**3

    vx_J[i+1] = vx_J[i] + ax*dt
    vy_J[i+1] = vy_J[i] + ay*dt
    x_J[i+1] = x_J[i] + vx_J[i+1]*dt
    y_J[i+1] = y_J[i] + vy_J[i+1]*dt
    
# Using Euler-Chromer method to simulate Earth's orbit due to Sun and Jupiter
for i in range(0, len(t)-1):
    # Calculate distance from Earth to Sun and from Earth to Jupiter
    r_E_to_Sun = np.sqrt(x_E[i]**2 + y_E[i]**2)
    r_E_to_J = myf.earth_jupiter_distance(x_E[i], y_E[i], x_J[i], y_J[i])
    
    # Calculate acceleration
    ax = -(G*M_sol* x_E[i]/r_E_to_Sun**3) + (G*M_J* (x_J[i]-x_E[i])/r_E_to_J**3)
    ay = -(G*M_sol* y_E[i]/r_E_to_Sun**3) + (G*M_J* (y_J[i]-y_E[i])/r_E_to_J**3)
    
    vx_E[i+1] = vx_E[i] + ax*dt
    vy_E[i+1] = vy_E[i] + ay*dt
    x_E[i+1] = x_E[i] + vx_E[i+1]*dt
    y_E[i+1] = y_E[i] + vy_E[i+1]*dt
    
    
fig = plt.figure(figsize=(8, 8))
plt.plot(0,0,'k.', label='Sun')
plt.plot(x_E, y_E, label="Earth's orbit")
plt.plot(x_J, y_J, label="Jupiter's orbit")
plt.xlabel('x-Position [AU]', fontsize=16)
plt.ylabel('y-Position [AU]', fontsize=16)
plt.title("Earth and Jupiter orbit over 10 Earth year", fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig('Q2a-Plot.png')

    
#%%
############################# Q2b ############################################

# Set Jupiter mass to be same as the mass of Sun
M_J = 1.0   # [M_sol]

# Build time array
dt = 0.0001             # [year]
end_time = 3.582        # [year]
t = np.arange(0,end_time,dt)

# Initialize zero arrays for Earth and Jupiter
x_E = np.zeros(len(t))
y_E = np.zeros(len(t))
vx_E = np.zeros(len(t))
vy_E = np.zeros(len(t))

x_J = np.zeros(len(t))
y_J = np.zeros(len(t))
vx_J = np.zeros(len(t))
vy_J = np.zeros(len(t))

# Initial Condition for Earth and Jupiter
x_E[0] = 1.0    # [AU]
y_E[0] = 0.0    # [AU]
vx_E[0] = 0.0   # [AU/yr]
vy_E[0] = 6.18  # [AU/yr]
x_J[0] = 5.2    # [AU]
y_J[0] = 0.0    # [AU]
vx_J[0] = 0.0   # [AU/yr]
vy_J[0] = 2.63  # [AU/yr]

# Using Euler-Cromer method to simulate Jupiter's orbit due to Sun 
# (assuming gravitational force due to Earth is small) first
for i in range(0, len(t)-1):
    r_J = np.sqrt(x_J[i]**2 + y_J[i]**2)
    
    # Calculate acceleration
    ax = -G * M_sol * x_J[i]/r_J**3
    ay = -G * M_sol * y_J[i]/r_J**3

    vx_J[i+1] = vx_J[i] + ax*dt
    vy_J[i+1] = vy_J[i] + ay*dt
    x_J[i+1] = x_J[i] + vx_J[i+1]*dt
    y_J[i+1] = y_J[i] + vy_J[i+1]*dt
    
# Using Euler-Chromer method to simulate Earth's orbit due to Sun and Jupiter
for i in range(0, len(t)-1):
    # Calculate distance from Earth to Sun and from Earth to Jupiter
    r_E_to_Sun = np.sqrt(x_E[i]**2 + y_E[i]**2)
    r_E_to_J = myf.earth_jupiter_distance(x_E[i], y_E[i], x_J[i], y_J[i])
    
    # Calculate acceleration
    ax = -(G*M_sol* x_E[i]/r_E_to_Sun**3) + (G*M_J* (x_J[i]-x_E[i])/r_E_to_J**3)
    ay = -(G*M_sol* y_E[i]/r_E_to_Sun**3) + (G*M_J* (y_J[i]-y_E[i])/r_E_to_J**3)
    
    vx_E[i+1] = vx_E[i] + ax*dt
    vy_E[i+1] = vy_E[i] + ay*dt
    x_E[i+1] = x_E[i] + vx_E[i+1]*dt
    y_E[i+1] = y_E[i] + vy_E[i+1]*dt
    
    
fig = plt.figure(figsize=(8, 8))
plt.plot(0,0,'k.', label='Sun')
plt.plot(x_E, y_E, label="Earth's orbit")
plt.plot(x_J, y_J, label="Jupiter's orbit")
plt.xlabel('x-Position [AU]', fontsize=16)
plt.ylabel('y-Position [AU]', fontsize=16)
plt.title("Earth and Jupiter (1000x mass) orbit", fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig('Q2b-Plot.png')


#%%
################################ Q2c #########################################

# Set Jupiter's mass back to normal 
M_J = 10e-3

# Build time array
dt = 0.0001             # [year]
end_time = 20           # [year]
t = np.arange(0,end_time,dt)

# Initialize zero arrays for asteroid and Jupiter
x_a = np.zeros(len(t))
y_a = np.zeros(len(t))
vx_a = np.zeros(len(t))
vy_a = np.zeros(len(t))

x_J = np.zeros(len(t))
y_J = np.zeros(len(t))
vx_J = np.zeros(len(t))
vy_J = np.zeros(len(t))

# Initial Condition for Asteroid and Jupiter
x_a[0] = 3.3    # [AU]
y_a[0] = 0.0    # [AU]
vx_a[0] = 0.0   # [AU/yr]
vy_a[0] = 3.46  # [AU/yr]
x_J[0] = 5.2    # [AU]
y_J[0] = 0.0    # [AU]
vx_J[0] = 0.0   # [AU/yr]
vy_J[0] = 2.63  # [AU/yr]

# Using Euler-Chromer method to simulate Jupiter's orbit due to Sun 
# (assuming gravitational force due to asteriod is small) first
for i in range(0, len(t)-1):
    r_J = np.sqrt(x_J[i]**2 + y_J[i]**2)
    
    # Calculate acceleration
    ax = -G * M_sol * x_J[i]/r_J**3
    ay = -G * M_sol * y_J[i]/r_J**3

    vx_J[i+1] = vx_J[i] + ax*dt
    vy_J[i+1] = vy_J[i] + ay*dt
    x_J[i+1] = x_J[i] + vx_J[i+1]*dt
    y_J[i+1] = y_J[i] + vy_J[i+1]*dt
    
# Using Euler-Cromer method to simulate Earth's orbit due to Sun and Jupiter
for i in range(0, len(t)-1):
    # Calculate distance from Earth to Sun and from Earth to Jupiter
    r_a_to_Sun = np.sqrt(x_a[i]**2 + y_a[i]**2)
    r_a_to_J = myf.earth_jupiter_distance(x_a[i], y_a[i], x_J[i], y_J[i])
    
    # Calculate acceleration
    ax = -(G*M_sol* x_a[i]/r_a_to_Sun**3) + (G*M_J* (x_J[i]-x_a[i])/r_a_to_J**3)
    ay = -(G*M_sol* y_a[i]/r_a_to_Sun**3) + (G*M_J* (y_J[i]-y_a[i])/r_a_to_J**3)
    
    vx_a[i+1] = vx_a[i] + ax*dt
    vy_a[i+1] = vy_a[i] + ay*dt
    x_a[i+1] = x_a[i] + vx_a[i+1]*dt
    y_a[i+1] = y_a[i] + vy_a[i+1]*dt

fig = plt.figure(figsize=(8, 8))
plt.plot(0,0,'k.', label='Sun')
plt.plot(x_a, y_a, label="Asteroid's orbit")
plt.plot(x_J, y_J, label="Jupiter's orbit")
plt.xlabel('x-Position [AU]', fontsize=16)
plt.ylabel('y-Position [AU]', fontsize=16)
plt.title("Asteroid and Jupiter orbit over 20 Earth year", fontsize=16)
plt.grid('on')
plt.legend(loc='lower right')
<<<<<<< HEAD
plt.savefig('Q2c-Plot.png')

=======
>>>>>>> e101bd8587053789d2321dbd8222c1181f7fd443
plt.show()
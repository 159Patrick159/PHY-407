# Contributors 
# Patrick Sandoval
# Kelvin Leong

################################# HEADER #####################################
# This python script holds the code, pseudo-code and plots for Q1 of Lab07 for 
# calculating trajectory of space garbage with adaptive and non-adaptive step
# Part of this code adapted from Nico Grisouard, Univ. of Toronto solution to 
# Newman 8.8, Space garbage.
##############################################################################
# Pseudocode:
# Import the necessary libraries
# Define the (force) function used in the RK4 method
# Calculate the x, vx, y, vy using the RK4 method with N = 1000 points
# To implement the adaptive method now
# Set initial h = 0.01, and target error-per-second delta; start the timer
# Initialize arrays to bookkeep x, vx, y, vy, t, and h points
# while t is less than end-time b:
    # Append all previous values to the bookkepping array, intialize rho
    # while rho < 1:
        # Perform the 2 times small step of RK4 using the current h (as x1)
        # Perform the 1 time big step of RK4 using the current h (as x2)
        # Calculate rho based on eq.4 in lab handhout
        # Calculate the new h using eq.8.52 in the textbook
        # If rho < 1, go back to beginning of this while-loop
    # Using the new h to perfrom 2 times small steps of RK4
# End the timer and calculate the time difference
# Plot the xy-trajectory from adaptive method and just the RK4 method
# Q1b:
# Calculate the trajectory again using just the RK4 method now with N = 10000
# points and time it (in order to achieve the same target error-per-second in
# adaptive method)
# Prints out the time it takes for both methods
# Q1c:
# Plot the bookkeeping h array vs t array

import numpy as np
import matplotlib.pyplot as plt
from time import time
from MyFunctions import rhs, RK4, adaptive_RK4


################################# Q1a ########################################
a = 0.0
b = 10.0
N = 1000  
h = (b-a)/N
r_rk4 = np.array([1., 0., 0., 1.], float)

# Results from non-adaptive method
xpoints_rk4, vxpoints_rk4, ypoints_rk4, vypoints_rk4, tpoints_rk4 = RK4(a,b,h,r_rk4)

# Adpative step
delta = 1e-6    # [m/s]
h = 0.01
r_adaptive = np.array([1., 0., 0., 1.], float)

# Start timing for adaptive method
start = time()
# Results from adaptive method
x_points_adaptive, vx_points_adaptive, y_points_adaptive, vy_points_adaptive, \
            tpoints_adaptive, h_adaptive = adaptive_RK4(a,b,h,delta,r_adaptive)
# End timing
end = time()
adaptive_time = end - start

# Plots xy-trajectory
plt.figure(figsize=(10,10))
plt.plot(xpoints_rk4, ypoints_rk4, 'k:', markersize=1.5, label='Non-adaptive')
plt.plot(x_points_adaptive, y_points_adaptive, 'r.', markersize=3.5, label='Adaptive')
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title('Trajectory of a ball bearing around a space rod.', fontsize=16)
plt.axis('equal')
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig('Q1a.png', dpi=150)

################################# Q1b ########################################
a = 0.0
b = 10.0
N = 10000  # let's leave it at that for now
h = (b-a)/N
r_rk4 = np.array([1., 0., 0., 1.], float)

start = time()
_, _, _, _, _ = RK4(a,b,h,r_rk4)
end = time()
non_adaptive_time = end - start

print("Time takes to achieve target error-per-second 10^-6 m/s for each method:")
print(f"Non-adaptive method: {non_adaptive_time: .5f} seconds")
print(f"Adaptive method: {adaptive_time: .5f} seconds")

################################ Q1c #########################################
# Plots the time step sizes as a function of time in adaptive method
plt.figure(figsize=(8,8))
plt.plot(tpoints_adaptive, h_adaptive, '.')
plt.xlabel("$t$ [s]")
plt.ylabel("Time steps $h$ in adaptive method [s]")
plt.grid()
plt.title("Size of the time step as a function of time for adaptive method", fontsize=14)
plt.tight_layout()
plt.savefig('Q1c.png')

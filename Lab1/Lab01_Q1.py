
# Collaborators:
# Patrick Sandoval
# Kelvin Leong

########################### HEADER ###############################
#  This python script contains the pseudo-code and explantion for Q1b the code, 
# plots and explanations for Q1c and code, plots and explanation for Q1d. 

############################################ Q1b ################################################

# Define/import needed physical constants from the scipy.constant library 
# Define the initial conditions for x0 and y0 and vx0  and vy0
# Define a time step (dt)
# Create and array for the time of the simulation with the prev. defined dt
# Create arrays of zeros for ax, ay, vx, vy, x, and y with the same length of the time domain
# Run a for loop from 1 to len(t) iterations
# Compute prev. instance r from x[i-1] and y[i-1]
# Compute the ith accelration from derived 1st order ode using the prev. position and distance instances
# Compute ith element of vx and vy using ith-1 vx and vy, the ith ax, ay, and dt
# Plot the positions x vs y and the vx vs t and vy vs time 

############################################ Q1c ################################################
# Import libraries
import numpy as np
import scipy.constants as c
import matplotlib.pyplot as plt


# Define initial conditions and physical constants
x0 = 0.47 # AU
y0 = 0.0 # AU
vx0 = 0.0 # AU/yr
vy0 = 8.17 # AU/yr
G = c.G # m3/s2/kg
G = G * (6.8459e-12)**3 * (3.154e7)**2 # AU3/yr2/kg
Ms = 1.98847e30 # kg

# Define timestep and time vector
dt = 0.0001
t = np.arange(0,1,dt)

# Define zero vectors for x,y,vx,vy,ax,ay
x = np.zeros(len(t))
y = np.zeros(len(t))
vx = np.zeros(len(t))
vy = np.zeros(len(t))
ax = np.zeros(len(t))
ay = np.zeros(len(t))

# Define initial conditions
x[0] = x0
y[0] = y0
vx[0] = vx0
vy[0] = vy0


# Run for loop from 1 to len(t) and compute ith elements
for i in range(1,len(t)):
    # Compute r from previous iteration
    r = np.sqrt(x[i-1]**2 + y[i-1]**2)
    # Compute acceleration
    ax[i] = -G*Ms*x[i-1]/r**3
    ay[i] = -G*Ms*y[i-1]/r**3
    # Compute velocity
    vx[i] = vx[i-1] + dt*ax[i]
    vy[i] = vy[i-1] + dt*ay[i]
    # Compute position
    x[i] = x[i-1] + dt*vx[i]
    y[i] = y[i-1] + dt*vy[i]


# Plot results
fig, (a0,a1,a2) = plt.subplots(figsize=(15,5),ncols=3)

a0.plot(x,y,c='k')
a0.set_title("Mercury's orbit",fontsize=18)
a0.set_xlabel("X-Position [AU]",fontsize=16)
a0.set_ylabel("Y-Position [AU]",fontsize=16)
a0.grid(ls='--')


a1.plot(t,vx,c='r')
a1.set_title("X-velocity component evolution",fontsize=18)
a1.set_ylabel(r"$V_x$",fontsize=16)
a1.set_xlabel("Time [yr]",fontsize=16)
a1.grid(ls='--')

a2.plot(t,vy,c='gray')
a2.set_title("Y-velocity component evolution",fontsize=18)
a2.set_ylabel(r"$V_y$",fontsize=16)
a2.set_xlabel("Time [yr]",fontsize=16)
a2.grid(ls='--')

plt.tight_layout()
plt.savefig("Q1c-Plot.png")


############################################ Q1d ################################################
# The code is exactly the same with the only difference being the (1+alpha/r^2) term in the 
# computation for the acceleration

alpha = 0.01 # AU2

# Run for loop from 1 to len(t) and compute ith elements
for i in range(1,len(t)):
    # Compute r from previous iteration
    r = np.sqrt(x[i-1]**2 + y[i-1]**2)
    # Compute acceleration
    ax[i] = -G*Ms*(1+alpha/r**2)*x[i-1]/r**3
    ay[i] = -G*Ms*(1+alpha/r**2)*y[i-1]/r**3
    # Compute velocity
    vx[i] = vx[i-1] + dt*ax[i]
    vy[i] = vy[i-1] + dt*ay[i]
    # Compute position
    x[i] = x[i-1] + dt*vx[i]
    y[i] = y[i-1] + dt*vy[i]

# Plot results
fig, (a0,a1,a2) = plt.subplots(figsize=(15,5),ncols=3)

a0.plot(x,y,c='k')
a0.set_title("Mercury's orbit",fontsize=18)
a0.set_xlabel("X-Position [AU]",fontsize=16)
a0.set_ylabel("Y-Position [AU]",fontsize=16)
a0.grid(ls='--')


a1.plot(t,vx,c='r')
a1.set_title("X-velocity component evolution",fontsize=18)
a1.set_ylabel(r"$V_x$",fontsize=16)
a1.set_xlabel("Time [yr]",fontsize=16)
a1.grid(ls='--')

a2.plot(t,vy,c='gray')
a2.set_title("Y-velocity component evolution",fontsize=16)
a2.set_ylabel(r"$V_y$",fontsize=14)
a2.set_xlabel("Time [yr]",fontsize=14)
a2.grid(ls='--')

plt.tight_layout()
plt.savefig("Q1d-Plot.png")
plt.show()






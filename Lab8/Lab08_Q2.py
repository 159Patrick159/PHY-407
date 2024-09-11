# Contributors 
# Patrick Sandoval
# Kelvin Leong

################################# HEADER ##########################################
# This python script holds the code, pseudo-code and plots for Q2 of Lab08 for    #  
# the simulation of shallow water using a flux conservative partial differetial   #
# equation of the 1D shallow water equation.                                      #
###################################################################################
# Pseudo-code
# Import needed libraries
# Define simulation parameters and physical constants
# Initialize solution arrays for mu and eta given grid spacing
# Define boundary conditions
# Define initial condition for mu and eta
# Start iterating through time with while loop
# Iterate through each particle spatially with for loop
# Use forward difference for lower bound
# Use backward difference for upper bound
# Use FTCS for everything else
# Overwrite prev. sol with current solution
# Increase time-step by dt
# Plot current solution for animation
# Check time for t=1,4 for saving fig
####################################### Q2b #######################################
import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import Feta, Fmu
from pylab import clf, plot, xlim, ylim, show, pause, draw

# Define boundary
L = 1 # m
J = 50
# Spatial spacing
del_x = L/J # m
# Temporal spacing 
del_t = 0.01 # s
# Constants
g = 9.8 # m/s^2
eta_b = 0
H = 0.01 # m
epsilon = del_t/1000

# Define solution arrays
mu = np.empty(J+1,float)
mup = np.empty(J+1,float)
eta = np.empty(J+1,float)
etap = np.empty(J+1,float)
# Boundary condtions
mu[0] = 0
mu[J] = 0
mup[0] = 0
mup[J] = 0
# Define times
t1, t2, t3 = 0, 1, 4 # s
# Define constants used for IC
A = 0.002 # m
mu0 = 0.5 # m
sigma = 0.05 # m

# Define initial conditions
x_domain = np.linspace(0,1,J+1)
mu[1:J] = 0
eta = H + A*np.exp(-(x_domain-mu0)**2/sigma**2) - np.mean(A*np.exp(-(x_domain-mu0)**2/sigma**2))

# Define constant 
C = del_t/2/del_x

# Define total time of solution
tend = 5
t = 0.0

plt.plot(x_domain,eta,marker='o',c='b',alpha=0.6)
plt.fill_between(x_domain,y1=0,y2=eta,color='blue',alpha=0.4)
plt.title(f"Shallow Water Simulation t=0s",fontsize=16)
plt.xlabel("Spatial Dimension (m)",fontsize=14)
plt.ylabel("Free Surface Altitude",fontsize=14)
plt.xlim([0,1])
plt.ylim([0,0.02])
plt.tight_layout()
plt.savefig("Q2_0s.pdf")


# Start time loop
while t < tend:
    # Start spatial loop
    for i in range(0,J+1):
        # Treat each boundary individually
        if i == 0:
            # Foward difference for lower bound
            mup[i] = mu[i] - C*(Fmu(mu[i+1],eta[i+1]) - Fmu(mu[i],eta[i]))
            etap[i] = eta[i] - C*(Feta(mu[i+1],eta[i+1]) - Feta(mu[i],eta[i]))
        elif i == J:
            # Backward difference for upper bound
            mup[i] = mu[i] - C*(Fmu(mu[i],eta[i]) - Fmu(mu[i-1],eta[i-1]))
            etap[i] = eta[i] - C*(Feta(mu[i],eta[i]) - Feta(mu[i-1],eta[i-1]))
        else:
            mup[i] = mu[i] - C*(Fmu(mu[i+1],eta[i+1])-Fmu(mu[i-1],eta[i-1]))
            etap[i] = eta[i] - C*(Feta(mu[i+1],eta[i+1])-Feta(mu[i-1],eta[i-1]))
    mu, mup = mup, mu
    eta, etap = etap, eta
    t += del_t
    clf()
    plt.title(f"Shallow Water Simulation t={round(t,3)}s",fontsize=16)
    plt.xlabel("Spatial Dimension (m)",fontsize=14)
    plt.ylabel("Free Surface Altitude",fontsize=14)
    plt.plot(x_domain,eta,label=r"$\eta$(x,t)",marker='o',c='b',alpha=0.6)
    plt.fill_between(x_domain,y1=0,y2=eta,color='b',alpha=0.4)
    plt.xlim([0,1])
    plt.ylim([0,0.02])
    plt.tight_layout()
    if abs(t-t2) < epsilon:
        plt.savefig("Q2_1s.pdf")
    if abs(t-t3) < epsilon:
        plt.savefig("Q2_4s.pdf")
    draw()
    pause(0.00001)






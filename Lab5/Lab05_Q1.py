# Contritbutors
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code, pseudo-code and plots
# for Q1 and all its subsections for Lab5
######################################## Q1.a ############################################
# Pseudocode 
# Import needed libraries
# Define physcal constants for the probelm and critical displacement
# Define the the intial conditions
# Define a time domain and a small time-step
# Call modified Euler_Cromer method for particle on spring
#   The only modification on this function is the different EOM
#   So it computes a different acceleration than the original function
# Plot results
##########################################################################################
# Import needed libraries
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import c
from MyFunctions import Euler_Cromer1D

# Define needed constants and critical displacement
m = 1 # kg
k = 12 # N/m
xc = c*np.sqrt(m/k)

# Define initial condition
x01 = 1 # m
x02 = xc # m
x03 = 10*xc
v0 = 0

# Define time array and time steps
dt = 0.001
t = np.arange(0,60,dt) # s

# Comppute the position of particles given inital positions
x1 = Euler_Cromer1D(t,x01,v0,m,k,dt)
x2 = Euler_Cromer1D(t,x02,v0,m,k,dt)
x3 = Euler_Cromer1D(t,x03,v0,m,k,dt)

fig, (a0,a1,a2) = plt.subplots(figsize=(8,6),nrows=3,\
                                sharex=True)

a0.plot(t,x1,c='k',label=r"$x_0$=1")
a0.set_title("Relativistic Particle on a Spring",fontsize=16)
a0.set_ylabel("X(t) (m)",fontsize=14)
a0.grid(ls='--')
a0.legend()

a1.plot(t,x2/1e7,c='r',label=r"$x_0$=$x_c$")
a1.set_ylabel(r"X(t) (10$^7$m)",fontsize=14)
a1.grid(ls='--')
a1.legend()

a2.plot(t,x3/1e8,c='gray',label=r"$x_0$=10$x_c$")
a2.set_ylabel("X(t) (10$^8$m)",fontsize=14)
a2.set_xlabel("Time [s]", fontsize=14)
a2.grid(ls='--')
a2.legend()

plt.subplots_adjust(hspace=0)
plt.savefig("Q1aplot.pdf")

######################################## Q1.b ############################################
# Pseudocode
# Take the FFT of the position arrays and shift the arrays such that zeroth freq is centered
# Normalize the results of the FFT by dividing by their max
# Generate the frequency domain from the time domain using np.fft.fftfreq
# Shift the frequecuency domain such that the zeroth frequency is centered
# Plot the results
##########################################################################################

# Start by taking the FFT of each signal
# Shift signal such that the zeroth frequncy is centered
G1 = np.fft.fftshift(np.fft.fft(x1))
G2 = np.fft.fftshift(np.fft.fft(x2))
G3 = np.fft.fftshift(np.fft.fft(x3))

# Normalize the fourier tranforms
G1 = np.abs(G1)/np.max(np.abs(G1))
G2 = np.abs(G2)/np.max(np.abs(G2))
G3 = np.abs(G3)/np.max(np.abs(G3))

# # Generate frequency domain from time domain and timestep
f = np.fft.fftshift(np.fft.fftfreq(len(t),dt))

# Check the periods of each wave
f1 = f[np.where(G1==np.max(G1))]

# Plot results
fig, a0 = plt.subplots(figsize=(8,4),ncols=1)
a0.plot(f*2*np.pi,G1,c='k',ls='--',label=r"$x_0$=1")
a0.plot(f*2*np.pi,G2,c='r',alpha=0.7,label=r"$x_0$=$x_c$")
a0.plot(f*2*np.pi,G3,c='k',alpha=0.8,label=r"$x_0$=10$x_c$")
a0.set_xlabel(r"Angular Frequency [rad/s]",fontsize=14)
a0.set_ylabel("Normalized Complex Amplitude",fontsize=14)
a0.set_title("FFT of Relativistic Particle on a Spring",fontsize=16)
a0.grid(ls='--')
a0.set_xlim([-5.,5.])
a0.legend(loc='upper left')
plt.tight_layout()
plt.savefig("Q1bplot.pdf")

######################################## Q1.c ############################################
# Pseudocode
# Find the harmonic frequecies of each signal
# by using np.where on the frequency arrays
# Take the inverse of these frequencies to find T
# Import the periods from the Gaussian quadrature
# method.
# Compute reisdual and plot Gaussian quadrature
# frequencies in FFT plot and compare results
##########################################################################################
# Compute periods from fourier 
Tf1 = 1/f[np.where(G1 == np.max(np.abs(G1)))[-1]]
Tf2 = 1/f[np.where(G2 == np.max(np.abs(G2)))[-1]]
Tf3 = 1/f[np.where(G3 == np.max(np.abs(G3)))[-1]]

# Initizlize periods from Lab03
Tclassical = 1.8102536 # s
T10xc = 11.6630122 # s
Txc = 2.12705948 #s

# Compute corresponding frequencies
fclassical = 1/Tclassical
f10xc = 1/T10xc
fxc = 1/Txc

# Plot results
fig, a0 = plt.subplots(figsize=(8,4),ncols=1)
a0.plot(f,G1,c='k',ls='--',label=r"$x_0$=1")
a0.plot(f,G2,c='r',alpha=0.7,label=r"$x_0$=$x_c$")
a0.plot(f,G3,c='k',alpha=0.8,label=r"$x_0$=10$x_c$")
a0.axvline(x=fclassical,ls='-.',c='orange',label='x=1')
a0.axvline(x=fxc,ls='-.',c='cyan',label=r'x=x$_c$')
a0.axvline(x=f10xc,ls='-.',c='gray',label=r'x=10x$_c$')
a0.set_xlabel(r"Linear Frequency [1/s]",fontsize=14)
a0.set_ylabel("Normalized Complex Amplitude",fontsize=14)
a0.set_title("FFT of Relativistic Particle on a Spring",fontsize=16)
a0.grid(ls='--')
a0.set_xlim([-0.7,0.7])
a0.legend()
plt.tight_layout()
plt.savefig("Q1cplot.pdf")


# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code and pseudo-code for Q2 and all its subsections 
# as well as the plots for the relativistic particle in a spring.

######################################## Q2.a ############################################
# Pseudo-code
# 
# Import needed libraries and functions to perform Gaussian quadrature integration
# Define needed physical quantities
# Define different values for N
# Define initial displacement
# Define boundaries
# Define variables to store sums
# Call gausswsab on the respective N's with the boundaries
# Initialize a for loop for one of the N's
#   Start summing up the weighted values of functions at the x's
#   provided by gausswxab
# Initialize another for loop for the remaining N
#   Sum up the weighted values of the functions at the x's
# Compute the theoretical classical period
# Compute fractional error for each result
# Print results

##########################################################################################
# Import needed libraries
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator

# Import gaussian quad functions
from MyFunctions import gaussxw, gaussxwab, v, Q2f

#Define needed constants
m = 1 # kg
k = 12 # N/m

# Define different values for N for integration
N1, N2 = 8, 16

# Define initial displacement for classical behaviour in (SI)
x0 = 0.01 #m

# Define boundaries for integration
a = 0
b = x0

# Compute the sample points and weights using gaussxwab
x1, w1 = gaussxwab(N1,a,b)
x2, w2 = gaussxwab(N2,a,b)
 
# Initialize integrals to 0
I1 = 0
I2 = 0

for i in range(N1):
    I1 += w1[i]*Q2f(x1[i],x0)
for j in range(N2):
    I2 += w2[j]*Q2f(x2[j],x0)

# Compute period's true value 
T_true = 2*np.pi*np.sqrt(m/k)
# Compute the fractional error
frac1 = abs(I1 - T_true)/I1
frac2 = abs(I2 - T_true)/I2


print("Classical Period:",T_true)
print("Gaussian Quad N=8,",I1,'Frac. Error:',frac1)
print("Gaussian Quad N=16,",I2,'Frac. Error:',frac2)
print()

######################################## Q2.b ############################################
# Pseudo-code
# 
# Initialize arrays to store the values of integrands and weigted integrands
# Set a for loop for one of the N's
#   Compute the integrand and weighted integrand for that N and store in array
# Set another for loop for the remaining N
#   Compute integrand and weighted integrand and stor in respective array
# Plot results

##########################################################################################
# Initialize arrays to populate with for loop
N1_int = np.zeros(N1)
N1_wint = np.zeros(N1)
N2_int = np.zeros(N2)
N2_wint = np.zeros(N2)

# Compute integrand and weigthed integrands
for i in range(N1):
    N1_int[i] = 4/v(x1[i],x0)
    N1_wint[i] = 4*w1[i]/v(x1[i],x0)
for j in range(N2):
    N2_int[j] = 4/v(x2[j],x0)
    N2_wint[j] = 4*w2[j]/v(x2[j],x0)

# Set minor tick direction
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (a0,a1) = plt.subplots(figsize=(12,4),ncols=2)
a0.plot(np.arange(1,N1+1,1),N1_int,marker='o',c='r',label=r'N=8: 4/g$_i$')
a0.plot(np.arange(1,N2+1,1),N2_int,marker='o',c='k',label=r'N=16: 4/g$_i$')
a0.set_title("Period Integrand",fontsize=16)
a0.set_xlabel("Number of Samples N",fontsize=14)
a0.set_ylabel(r"Intgrand ( 4/g$_i$ )",fontsize=14)
a0.xaxis.set_minor_locator(MultipleLocator(0.2))
a0.yaxis.set_minor_locator(MultipleLocator(20))
a0.yaxis.set_ticks_position('both') 
a0.xaxis.set_ticks_position('both')
a0.legend()

a1.plot(np.arange(1,N1+1,1),N1_wint,marker='o',c='r',label=r'N=8: 4w$_i$/g$_i$')
a1.plot(np.arange(1,N2+1,1),N2_wint,marker='o',c='k',label=r'N=16: 4w$_i$/g$_i$')
a1.set_title("Period Weighted Integrand",fontsize=16)
a1.set_xlabel("Number of Samples N",fontsize=14)
a1.set_ylabel(r"Weigthed Intgrand ( 4w$_i$/g$_i$ )",fontsize=14)
a1.xaxis.set_minor_locator(MultipleLocator(0.2))
a1.yaxis.set_minor_locator(MultipleLocator(0.01))
a1.yaxis.set_ticks_position('both') 
a1.xaxis.set_ticks_position('both')
a1.legend()

plt.tight_layout()
plt.savefig("Q2bPlot.png")

######################################## Q2.c ############################################
# Please refer to the written report for the full derivation

######################################## Q2.d ############################################
# Pseudo-code
# 
# Define new number of samples
# Get weights and x's from gausswx
# Define summation variable
# Start for loop
#   Compute weighted sum 
# Compute fractional error
# Print results

##########################################################################################

# We will implement the same code from Q2a for a number of samples of N=200 for a small 
# amplitude relativistic harmonic oscillator.

# Define number of samples
N3 = 200 

# Use gaussxwab to get weights and x values for integration
x3, w3 = gaussxwab(N3,a,b)

# Sum over weighted elements 
I3 = 0
for i in range(N3):
    I3 += w3[i]*Q2f(x3[i],x0)

# Compute fractional error
frac3 = abs(I3-T_true)/I3
# Print results
print("For 200 Samples")
print("Classical Period of Oscillator:",T_true)
print("Gaussian Quadrature:",I3,"Frac. Error:",frac3)
print()

######################################## Q2.e ############################################
# Pseudo-code
# 
# Compute the critical displacement derived from prev. question
# Initialze an array for x0's
# Define the number of samples for Gaussina quad procedure
# Call gausswa for the specified numbers
# Define empty array to store periods
# Iterate through x0
# Set temporary variable for sum
# Define new x and w for x0_ith upper bound
# Perform weighted sum
# Append result to period array
# Plot results

##########################################################################################
# Define speed of light
c = 2.998e8 # m/s

# Define critical displacement
xc = c*np.sqrt(m/k)

# Create array of x0s
x0a = np.linspace(1,10*xc,500)

# Define number of elements
N = 500

# Define x and w for N
x, w = gaussxw(N)

# Create array for periods
T = []

# Compute T for each x0
for val in x0a:
    # Set temp period
    T_tmp = 0

    # Adjust x and w for new bounds
    x4 = 0.5*(val)*x+0.5*(val)
    w4 = 0.5*(val)*w

    for i in range(N):
        T_tmp += w4[i]*Q2f(x4[i],val)
    T.append(T_tmp)

# Plot results
fig, a0 = plt.subplots(figsize=(8,5),ncols=1)
a0.plot(x0a,T,c='k')
a0.set_title("Period of Relativistic Harmonic Oscillator for Large Applitudes",fontsize=16)
a0.set_xlabel(r"Amplitude x$_0$ [m]",fontsize=14)
a0.set_ylabel("Period T [s]",fontsize=14)
a0.xaxis.set_minor_locator(MultipleLocator(0.1e8))
a0.yaxis.set_minor_locator(MultipleLocator(0.2))
a0.yaxis.set_ticks_position('both') 
a0.xaxis.set_ticks_position('both')
plt.tight_layout()
plt.savefig("Q2e.png")
plt.show()

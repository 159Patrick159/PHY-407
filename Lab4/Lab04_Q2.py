# Contributors:
# Patrick Sandoval
# Kelvin Leong

################################ HEADER ######################################
# This file contains the main body code, pseudocode, and print out for Q2
# and all its consequent sub-questions.
##############################################################################

############################### Q2b ##########################################
# Import needed libraries
from xml.dom.xmlbuilder import DOMInputSource
import numpy as np
from numpy.linalg import eigvalsh
from scipy.integrate import simps
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from MyFunctions import make_wf

# Define function for H's matrix element
def Hmatrix(m,n):
    from scipy.constants import hbar
    '''Retruns the corrsponding matrix element given the 
    labs definition of the Hamiltonian matrix'''
    # Define problem constants in given units
    a = 10 # eV
    M = 9.1094e-31 # kg
    L = 5 # AA

    # Convert all constants to SI
    L = L*1e-10 # m
    a = a*1.60218e-19 # Joules  

    # Check for intial condition if m = n
    if m == n:
        return(0.5*a + (np.pi*hbar*m)**2/2/M/L**2)
    # Check for second condition of m != n
    if (m != n):
        # Check for even m
        if (m%2 == 0):
            # m is even but n is odd
            if (n%2 != 0):
                return(-8*a*m*n/np.pi**2/(m**2-n**2)**2)
            # Both m n are even
            else:
                return(0)
        # Check for odd case
        else:
            # m is odd but n is even
            if (n%2 == 0):
                return(-8*a*m*n/np.pi**2/(m**2-n**2)**2)
            # Both m and n are odd
            else:
                return(0)

############################### Q2c ##########################################
# Deine size of matrix
mmax = 10
nmax = 10

# Define matrix with specified dimensions 
H=np.zeros((mmax,nmax))

# Compute the matrix elements
for m in range(1,mmax+1):
    for n in range(1,nmax+1):
        H[m-1,n-1] = Hmatrix(m,n)

# Calculate the eigenvalues of the Hamiltonian
Eig = eigvalsh(H) # in Joules

# Convert from Joules to eV
Eig = Eig*6.242e+18 # eV
print(f"Hamiltonian ({mmax}x{nmax}) Matrix")
print('='*70)
print("Energy Eigenvalues")
print(Eig)
print()

############################### Q2d ##########################################
# Define size of hamiltonian
mmax2, nmax2 = 100,100
# Initialize empty hamiltonian
H100 = np.zeros((mmax2,nmax2))

# Compute the hamiltonian
for m in range(1,mmax2+1):
    for n in range(1,nmax2+1):
        H100[m-1,n-1] = Hmatrix(m,n)

# Calculate eigenvalues
Eig100 = eigvalsh(H100)

# Convert to eV
Eig100 = Eig100*6.242e+18 # eV

print(f"Hamiltonian ({mmax2}x{nmax2}) Matrix")
print('='*70)
print("Energy Eigenvalues")
print(Eig100[:10])
print()

############################### Q2e ##########################################
# Pseuod-code
# Find the eigenvalues and eigenvectors of the (100x100) Hamiltonian
# From numpy documentation the ith columns correspond to the ith 
# eigenvalues
# Transpose the eigenvector matrix
# Select 1st, 2nd and 3rd row for ground, 1st and 2nd excited states
# Generate the spatial domain for the length of the well
# Define function thats supeimposes the elements of the eigenvector
# to compute the wavefunction at that energy state
# Call make_wf on corresponding eigenvector and generated domain
# Integrate the square norm of wave function
# Divide by sqrt of integration results to normalize wf
# Plot the probability amplitude of the wf
##############################################################################
# Find the eigenvalue and eigenvector for the (100x100) Hamiltonian
# Each column of V is an eigenvector of H
x, V = np.linalg.eigh(H100) # SI units
# Transpose matrix so rows are the eigenvectors
V = V.T
# Select the first three eigenvectors corresponding to ground, first and second excited state
V0, V1, V2 =  V[0], V[1], V[2]
x0, x1, x2 = x[0], x[1], x[2]

# Generate domain
LA = 5 # AA
domain = np.linspace(0,LA,1000) # AA
# Compute wavefunction
WF0 = make_wf(V0,domain)
WF1 = make_wf(V1,domain)
WF2 = make_wf(V2,domain)

# Normalize the wavefunctions
N0 = simps(np.power(np.abs(WF0),2))
N1 = simps(np.power(np.abs(WF1),2))
N2 = simps(np.power(np.abs(WF2),2))

WF0 = WF0/np.sqrt(N0)
WF1 = WF1/np.sqrt(N1)
WF2 = WF2/np.sqrt(N2)

# Set minor tick direction
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

# Plot probability amplitudes
fig, a0 = plt.subplots(figsize=(8,6),ncols=1)
a0.plot(domain,WF0**2,c='k',label='Ground State')
a0.plot(domain,WF1**2,c='r',label=r'1$^{st}$ Excited State')
a0.plot(domain,WF2**2,c='gray',label=r'2$^{st}$ Excited State')
a0.set_xlabel(r"Spatial Dimension $\AA$",fontsize=16)
a0.set_ylabel(r"Probability Amplitude $|\Psi|^2$",fontsize=16)
a0.set_title("Probability Amplitude in Asymmetric Well",fontsize=18)
a0.xaxis.set_minor_locator(MultipleLocator(0.2))
a0.yaxis.set_minor_locator(MultipleLocator(0.0001))
a0.yaxis.set_ticks_position('both') 
a0.xaxis.set_ticks_position('both')
a0.legend()
plt.savefig("Q2.png")



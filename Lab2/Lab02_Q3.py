# Collaborators:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code and pseudocode for all the subquestion in Q2 of Lab02.
# Additionally we explain some parts of the code and provide useful comments. This code
# compares Simpson's rule and Scipy.integrate.quad rotuine by calculating the Stefan-Boltzmann
# constant and testing its relative error.

######################################## Q3.b ############################################
# Pseudo-code for integrating Plank's law.
# Define intrgrand as shown in eq.8 of lab manual
# Define a function that will calculate the spectral radiance
# for a given T, and upper and lower bounds
# 
# Import simpsons integrating function from MyFunctions.py
# Import needed constants
# Compute the constant C1 using the imported values
# Perform the integration using Simpson's rule
# Multiply C1 with integrated value and return that

######################################## CODE ########################################

# Import needed libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
from scipy.optimize import curve_fit
import scipy.constants as u

# Define function that will be integrated
def f(x):
    return(x**3/(np.exp(x)-1))


def BB_Radiance(x,T,upper,lower):
    '''This function computes the blackbody radiance
    for a given temperature by integrating the Planks law'''
    from MyFunctions import simpsons
    from scipy.integrate import quad
    import scipy.constants as u
    # Define needed constants
    k = u.k
    h = u.h
    c = u.c

    # Define C1
    C1 = (2*np.pi*k**4*T**4)/(h**3*c**2)

    # Because the lim of f(x) as x-> inf = 0
    # Then setting a large upper bound would 
    # be a good approximation.
    
    # Compute integral
    integral1 = simpsons(f,lower,upper,len(x))
    integral2 = quad(f,lower,upper)[0]

    # Return radiance
    return(C1*integral1,C1*integral2)

######################################## Q3.c ############################################   
# Generate domain for x with upper and lower bounds
up = 200
lo = 1e-5
N = 50
xdomain = np.linspace(lo,up,N)

# Generate Temperature domain
Tdomain = np.linspace(10,10000,N)

# Compute radiance for every T
Ws = np.zeros(N)
Wq = np.zeros(N)
for i,T in enumerate(Tdomain):
    Ws[i] = BB_Radiance(xdomain,T,up,lo)[0]
    Wq[i] = BB_Radiance(xdomain,T,up,lo)[1]

# Compute T^4
T4 = Tdomain**4

# Slope of W/T4 is Steffan-Boltzmann constant
def linear(x,m,c):
    return(x*m + c)

popts,pcovs = curve_fit(linear,T4,Ws,absolute_sigma=True)
poptq,pcovq = curve_fit(linear,T4,Wq,absolute_sigma=True)

# Print results
print("Stefan-Boltzmann Constant from Simpson Routine:",popts[0])
print("Stefan-Boltzmann Constant from scipy.integrate.quad Routine:",poptq[0])
print("Stefan-Boltzmann Constant true value:",u.sigma)
print()
# Compute relative error
rels = (popts[0] - u.sigma)/u.sigma
relq = (poptq[0] - u.sigma)/u.sigma
print("Simpson relative error:",rels)
print("Quad relative error:",relq)

# Set plotting style
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig,a0 = plt.subplots(figsize=(8,5),ncols=1)

a0.scatter(T4,Ws,c='k',label="Simpson Rule",zorder=4)
a0.scatter(T4,Wq,c='k',alpha=0.5,marker="^",label="Quadpack",zorder=3)
a0.plot(T4,linear(T4,*popts),c='r',label='Least-square fit',ls='--',zorder=2)
a0.plot(T4,linear(T4,*poptq),c='r',ls='--',zorder=1)

a0.set_xlabel(r"Quartic Temperature [K$^4$]",fontsize=14)
a0.set_ylabel(r"Radiant Emittance [Watt $\cdot$ m$^{-2}$]",fontsize=14)
a0.set_title("Estimating Stefan-Boltzmann Constant",fontsize=16)
a0.xaxis.set_minor_locator(MultipleLocator(0.02e16))
a0.yaxis.set_minor_locator(MultipleLocator(0.1e8))
a0.legend()
a0.yaxis.set_ticks_position('both') 
a0.xaxis.set_ticks_position('both')
plt.tight_layout()
plt.savefig("Q3cPlot.png")
plt.show()





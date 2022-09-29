# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code and pseudo-code for Q2 and all its subsections 
# as well as the plots for the relativistic particle in a spring.

######################################## Q2.a ############################################
# Import needed libraries
from random import gauss
import re
import numpy as np
# Import gaussian quad functions
from MyFunctions import gaussxw, gaussxwab, v

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

# Define function to integrate
def f(x):
    return(4/v(x))
 
# Initialize integrals to 0
I1 = 0
I2 = 0

for i in range(N1):
    I1 += w1[i]*f(x1[i])
for j in range(N2):
    I2 += w2[i]*f(x2[i])

print("Classical Period:",2*np.pi*np.sqrt(m/k))
print("Gaussian Quad N=8,",I1)
print("Gaussian Quad N=16,",I2)
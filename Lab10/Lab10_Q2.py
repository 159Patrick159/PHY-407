########################## 
# Contributors           # 
# Patrick Sandoval       #
# Kelvin Leong           #
##########################

################################# HEADER ##########################################
# This python script holds the code for question 2 of lab 10 for PHY407           #  
# we implement Montecarlo integration techniques to compute the volume of a       #
# 10-dimensional hyperspher                                                       #
###################################################################################
# Import libraries
from random import random
import numpy as np

# Define function to integrate
def f(r):
    if np.sum(r**2) <= 1:
        return(1)
    else:
        return(0)

# Define bounds of integral
a = -1
b = 1


# Define sample size
N = 1e6
# Define number of dimensions
dof = 10
# Define sum variable
sum = 0

# Define lenght of hypercube enclosing hypersphere
L = 2**dof

for i in range(int(N)):
    # Generate 10 randomom number from -1 to 1
    xi = np.array([(b-a)*random() + a for _ in range(dof)])
    # Compute function for random numbers
    y = f(xi)
    sum += y

result = L/N*sum
print(f"Volume of sphere in {dof} dimensions:",result)
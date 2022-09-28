# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code and pseudo-code for Q1 and all its subsections 
# as well as the plots and numerical values of our differentiation schemes.

######################################## Q1.a ############################################
# Pseudo-code:
# Import error function from MyFunction.py script
# Define array for h starting from 1e-16 to 1 in steps of 10
# /
# Import error function from MyFunctions.py
from MyFunctions import erf
import numpy as np
h = [10**val for val in range(-16,0)]
print(h)

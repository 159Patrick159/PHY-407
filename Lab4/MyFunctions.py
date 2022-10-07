# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the functions for Lab04 in the process of answering 
# questions Q2 and Q3 and their subquestions.
##########################################################################################
import numpy as np

def make_wf(EigV,x):
    '''Given the eigenvector it calculates
    the corresponding wavefunction by the superposition
    principle'''
    LA = 5 # AA
    tmp = []
    for i in range(1,len(EigV)-1):
        tmp.append(EigV[i-1]*np.sin(i*np.pi*x/LA))
    return(sum(tmp))

# Function for q3 a and b
def Q3ab_f(x,c):
    return 1 - np.exp(-c*x)

# First derivative of function for q3 a and b
def Q3ab_f_derivative(x,c):
    return c*np.exp(-c*x)

# Relaxation method that numerically search for solution of non-linear equations
def relaxation_method(f, fprime, x, c, accuracy):
    """
    Parameters: f : function
                fprime : 1st derivative of function
                x : Initial guess of solution
                c : A parameter for function f
                accuracy : Desired accuracy of the numerical solution
    Returns: x1 : The numerical solution with desired accuracy
             error : Error on the numerical solution
             counter : Number of steps takes to find x1
    """
    x1 = x
    error = 1.0
    counter = 0
    # Loop until error is smaller than target accuracy
    while error > accuracy:
        x1, x2 = f(x1,c), x1
        error = abs((x1 - x2)/(1- 1/fprime(x1,c)))
        #print(x1, error)
        counter += 1
    return x1, error, counter

# Overrelaxation method
def overrelaxation_method(f, fprime, x, c, accuracy, omega):
    """
    Parameters: f : function
                fprime : 1st derivative of function
                x : Initial guess of solution
                c : A parameter for function f
                accuracy : Desired accuracy of the numerical solution
                omega: Parameter for overrelaxation
    Returns: x1 : The numerical solution with desired accuracy
             error : Error on the numerical solution
             counter : Number of steps takes to find x1
    """
    x1 = x
    error = 1.0
    counter = 0
    # Loop until error is smaller than target accuracy
    while error > accuracy:
        x1, x2 = (1+omega)*f(x1,c) - omega*x1, x1
        error = abs((x1 - x2)/(1- 1/((1+omega)*fprime(x1,c) - omega)))
        #print(x1, error)
        counter += 1
    return x1, error, counter

# Function for q3 c
def Q3c_f(x):
    return 5*np.exp(-x) + x - 5

# First derivative of function for q3 c
def Q3c_f_derivative(x):
    return -5*np.exp(-x) + 1

# Binary search method (procedure in pg.264)
def binary_search(f, x1, x2, accuracy, binary_counter):
    """
    Parameters: f : function
                x1 : Initial guess of lower bound of solution
                x2 : Initial guess of upper bound of solution
                accuracy : Desired accuracy of the numerical solution
                binary_counter: Dummy variable used for return
    Returns: Recursion! For numerical solution 
             binary_counter : Number of steps takes to find solution / Depth 
                              of recursion
    """
    # Use recursion to implement binary search
    binary_counter += 1
    
    f_x1 = f(x1)
    #f_x2 = f(x2)
    
    # Calculate the midpoint
    x_half = 0.5 * (x1 + x2)
    f_x_half = f(x_half)
    
    # Set x1 to midpoint if f(x_mid) has same sign as f(x_1)
    # else set x2 to midpoint if f(x_mid) has same sign as f(x_2)
    if (f_x1 * f_x_half > 0):
        x1 = x_half
    else:
        x2 = x_half
        
    # Recursion step:
    # If |x1-x2| less than target accuracy, take the midpoint as solution
    # else calls binary_search()
    if(abs(x1 - x2) < accuracy):
        return 0.5*(x1 + x2), binary_counter    #Base case
    else:
        return binary_search(f, x1, x2, accuracy, binary_counter)
    
# Newton's method
def Newton_method(f, fprime, x, accuracy):
    """
    Parameters: f : function
                fprime : 1st derivative of function
                x : Initial guess of solution
                accuracy : Desired accuracy of the numerical solution
    Returns: x: The numerical solution with desired accuracy
             counter : Number of steps takes to find x
    """
    counter = 0
    delta = 1.0
    while abs(delta) > accuracy:
        delta = f(x)/fprime(x)
        x -= delta
        counter +=1
    return x, counter
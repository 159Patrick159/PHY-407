# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds all the functions for lab03 in the process of answering 
# all the questions and subquestions.
import numpy as np

# Define error funtion for Q1
def erf(x):
    return(np.exp(-np.power(x,2)))

# Define the analytic first derivative of error function
def first_derivative_erf(x):
    return -2*x*np.exp(-np.power(x,2))

# Numerical derivative of forward difference (Eq. 5.90)
def fwd_derivative(f, x, h):
    """
    Parameters: f: Function to take derivative on
                x: Point of interest to take the derivative on
                h: Step size
    Returns: The first derivative on x using forward difference method
    """
    return (f(x+h) - f(x))/h

# Numerical derivative of central difference (Eq. 5.98)
def central_derivative(f, x, h):
    """
    Parameters: f: Function to take derivative on
                x: Point of interest to take the derivative on
                h: Step size
    Returns: The first derivative on x using central difference method
    """
    return (f(x+0.5*h) - f(x-0.5*h))/h



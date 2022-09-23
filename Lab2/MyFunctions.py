#!/usr/bin/env python3
# -*- coding: utf-8 -*-

############################ HEADER ##########################################
# This python file contains our self-defined function used for Lab02
##############################################################################

import numpy as np

def find_mean(data):
    '''Returns the mean of a data set'''
    import numpy as np
    return(np.sum(data)/len(data))

# Define integrand in the integral in Q2 as a function
def q2_f(x):
    return 4/(1+x**2)
    
# Calculate area under the curve using trapezoidal rules (eq. 5.3)
def trapezoidal(f,a,b,N):
    """
    Parameters: f: [function(float): float] Integrand function
                a: Lower bound of the integral
                b: Upper bound of the integral
                N: Number of slices

    Returns: Trapezoidal approximation of the integral
    """
    # Width of slice
    h = (b - a) / N
    
    result = 0.5 * (f(a) + f(b))
    for k in range(1, N):
        result += f(a + k*h)
    result = h * result
    return result

# Calculate area under the curve using simpson's rules (eq. 5.9)
def simpsons(f,a,b,N):
    """
    Parameters: f: [function(float): float] Integrand function
                a: Lower bound of the integral
                b: Upper bound of the integral
                N: Number of slices

    Returns: Simpson's approximation of the integral
    """
    # Width of slice
    h = (b - a) / N
    
    odd_sum = 0
    even_sum = 0
    
    for k in range(1, N):
        if (k % 2)== 1:     # k is odd
            odd_sum += f(a + k*h)
        else:
            even_sum += f(a + k*h)
    result = h/3 * (f(a) + f(b) + 4 * odd_sum + 2 * even_sum)
    return result


# Question 4 functions
# Polynomial function in q4 (simplified form)
def p(u):
    return (1-u)**8

# Polynomial function in q4 (expanded form)
def q(u):
    return 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

# Calculate the roundoff error of p(u)-q(u)
def q4b_roundoff_error(u):
    # p(u)-q(u) can be written as a sum of 10 terms, each of these 10 terms
    # is power of (1-u) or of u
    term_p = ((1-u)**8)
    term_q0 = -1
    term_q1 = (8*u)
    term_q2 = (-28*u**2)
    term_q3 = (56*u**3)
    term_q4 = (-70*u**4)
    term_q5 = (+56*u**5)
    term_q6 = (-28*u**6)
    term_q7 = (+8*u**7)
    term_q8 = (-u**8)

    # Take the square of each of these term and sum them over
    sum_square_x = term_p**2 + term_q0**2 + term_q1**2 + term_q2**2 + term_q3**2 + term_q4**2 \
                    +term_q5**2 +term_q6**2 +term_q7**2 + term_q8**2
        
    # Calculate the root mean square of these 10 terms
    rms = np.sqrt(sum_square_x/10)
    
    C = 1e-16
    return C * np.sqrt(10) * rms

# Calculate the fractional error p(u)-q(u)
def q4c_fractional_error(u):
    # Same roundoff error calculation same as previous function
    term_p = ((1-u)**8)
    term_q0 = -1
    term_q1 = (8*u)
    term_q2 = (-28*u**2)
    term_q3 = (56*u**3)
    term_q4 = (-70*u**4)
    term_q5 = (+56*u**5)
    term_q6 = (-28*u**6)
    term_q7 = (+8*u**7)
    term_q8 = (-u**8)
    sum_square_x = term_p**2 + term_q0**2 + term_q1**2 + term_q2**2 + term_q3**2 + term_q4**2 \
                    +term_q5**2 +term_q6**2 +term_q7**2 + term_q8**2
    rms = np.sqrt(sum_square_x/10)
    C = 1e-16
    rdf_error = C * np.sqrt(10) * rms

    # Calculate the sum of all 10 terms
    sum_x = term_p+term_q0+term_q1+term_q2+term_q3+term_q4+term_q5+term_q6+term_q7+term_q8
    
    # Calculate the fractional error according to eq(4)
    fractional_error = rdf_error/sum_x
    return fractional_error

# Product and quotient function in q4d
def q4d_f(u):
    return u**8/((u**4) * (u**4))
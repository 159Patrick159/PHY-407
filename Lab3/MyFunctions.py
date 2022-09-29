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

# Define Gaussian quad functions from textbook 
def gaussxw(N):
    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

# Define velocity function for Q2
def v(x,x0):
    #Define needed constants
    m = 1 # kg
    k = 12 # N/m
    c = 2.998e8 # m/s
    # Compute term
    term = c*np.sqrt((k*(x0**2-x**2)*(2*m*c**2 + k*(x0**2-x**2)/2))/(2*(m*c**2 + k*(x0**2-x**2)/2)**2))
    return(term)

# Define function to integrate for Q2
def Q2f(x,x0):
    return(4/v(x,x0))
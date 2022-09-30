# Contributors:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds all the functions for lab03 in the process of answering 
# all the questions and subquestions.
##########################################################################################
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



# The Hermite Polynomials in Eq.9
def Hermite_Polynomial(n,x,print_last_two=False):
    """
    Parameters: n: nth energy level, integer greater equal than 0
                x: Position
                print_last_two: (Optional) If True, will return last two Hermite 
                                polynomials; Default=False
    Return: nth Hermite Polynomial
            (if print_last_two==True, nth and (n-1)th Hermite Polynomials)
    """
    # Define the base cases H0 and H1
    # If argument x is a float, H0=1; if argument x is an array, H0=Array of 1
    if type(x)==float or type(x)==np.float64:
        H_0 = 1
    else:
        H_0 = np.ones(len(x))
    H_1 = 2*x
    
    # Put the two base cases in the array
    H_array = [H_0, H_1]
    
    # Calculate next Hermite Polynomial from the previous and this polynomial
    # Then append to array
    for i in range(1, n):
        H_next = 2*x*H_array[i] - 2*i*H_array[i-1]
        H_array.append(H_next)
        
    # If print_last_two==True, return nth and (n-1)th Hermite Polynomials
    # Specifically used to calculate the first derivative of Hermite Polynomials
    if print_last_two:
        if n == 0: return H_0, 0    # if n=0 case, then return H0 and 0
        return H_array[-1], H_array[-2]
        
    # Return the nth Hermite Polynomial
    if n == 0: return H_0           # if n=0 case, then return H0
    return H_array[-1]


# Define the nth energy level wavefunction from Eq.8
def Wave_function(n, x):
    """
    Parameters: n: nth energy level, integer greater equal than 0
                x: Position
    Return: nth wavefunction
    """
    H_n = Hermite_Polynomial(n, x)
    return np.exp(-x**2/2)/np.sqrt(2**n * np.math.factorial(n) * np.sqrt(np.pi)) * H_n


# Define the nth energy level first derivative of wavefunction from Eq.11
def Wave_function_1st_derivative(n, x):
    """
    Parameters: n: nth energy level, integer greater equal than 0
                x: Position
    Return: 1st derivative of nth wavefunction
    """
    H_n, H_n_minus1 = Hermite_Polynomial(n, x, print_last_two=True)
    return np.exp(-x**2/2)/np.sqrt(2**n * np.math.factorial(n) * np.sqrt(np.pi)) \
        * (-x*H_n + 2*n*H_n_minus1)


# Integrand in Eq.12 but with change of variables
def integrand_expected_position_square(n,z):
    """
    Parameters: n : nth energy level, integer greater equal to zero
                z : Position domain after change of variable
    Return: Integrand of <x^2> after change of variable
    """
    Psi_n = Wave_function(n, np.tan(z))
    return (np.tan(z)/np.cos(z))**2 * (abs(Psi_n))**2

# Integrand in Eq.13 but with change of variables
def integrand_expected_momentum_square(n,z):
    """
    Parameters: n : nth energy level, integer greater equal to zero
                z : Position domain after change of variable
    Return: Integrand of <p^2> after change of variable
    """
    dPsi_n = Wave_function_1st_derivative(n, np.tan(z))
    return (abs(dPsi_n)/np.cos(z))**2
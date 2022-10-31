# Contributors 
# Patrick Sandoval
# Kelvin Leong

############################### HEADER #######################################
# This file contains all the self-written function used for Lab7
##############################################################################

# Impoer libraries
import numpy as np

# Constants
m = 9.1094e-31
hbar = 1.0546e-34
e = 1.6022e-19
epsilon_0 = 8.854e-12 # F/m
a = 5e-11 # m (Bohr radius)
h = 0.0001*a # step size
rinf = 20*a


# Potential function
def V(r):
    return -e**2/4/np.pi/epsilon_0/r

def f(v,r,E,l):
    R = v[0]
    S = v[1]
    fR = S/r**2
    fS = ((2*m*r**2/hbar**2)*(V(r)-E)+l*(l+1))*R
    return np.array([fR,fS] ,float)

# Calculate the wavefunction for a particular energy
def solve_En(E,l):
    psi = 0.0
    phi = 1.0
    v = np.array([psi,phi] ,float)
    for r in np.arange(h,rinf,h):
        k1 = h*f(v,r,E,l)
        k2 = h*f(v+0.5*k1,r+0.5*h,E,l)
        k3 = h*f(v+0.5*k2,r+0.5*h,E,l)
        k4 = h*f(v+k3,r+h,E,l)
        v += (k1+2*k2+2*k3+k4)/6
    return v[0]

def solve_Rn(E,l):
    R = []
    S = 1.0
    v = np.array([0.0,S] ,float)
    for r in np.arange(h,rinf,h):
        R.append(v[0])
        k1 = h*f(v,r,E,l)
        k2 = h*f(v+0.5*k1,r+0.5*h,E,l)
        k3 = h*f(v+0.5*k2,r+0.5*h,E,l)
        k4 = h*f(v+k3,r+h,E,l)
        v += (k1+2*k2+2*k3+k4)/6
    return(R)

def En(n):
    return(-15*e/n**2,-13*e/n**2)

# The right-hand-side of the equations for Q1
def rhs(r):
    """ 
    INPUT:
    r = [x, vx, y, vy] are floats (not arrays)
    note: no explicit dependence on time
    OUTPUT:
    1x2 numpy array, rhs[0] is for x, rhs[1] is for vx, etc"""
    M = 10.
    L = 2.

    x = r[0]
    vx = r[1]
    y = r[2]
    vy = r[3]

    r2 = x**2 + y**2
    Fx, Fy = - M * np.array([x, y], float) / (r2 * np.sqrt(r2 + .25*L**2))
    return np.array([vx, Fx, vy, Fy], float)

# 4th order Runge-Kutta method, adapted from Newman's odesim.py 
def RK4(a,b,h,r_init):
    """
    Parameters: a: Start time of simulation
                b: End time of simulation
                h: Time step (Constant for RK4)
                r_init: Initial conditions in the order of [x, vx, y, vy]

    Returns: xpoints_rk4 : Array that stores the x-position of particle
            vxpoints_rk4 : Array that stores the x-velocity of particle
             ypoints_rk4 : Array that stores the y-position of particle
            vypoints_rk4 : Array that stores the y-velocity of particle
             tpoints_rk4 : Array that stores time after each time step
    """
    tpoints_rk4 = np.arange(a, b, h)
    xpoints_rk4 = []
    vxpoints_rk4 = []  # the future dx/dt
    ypoints_rk4 = []
    vypoints_rk4 = []  # the future dy/dt
    
    # below: ordering is x, dx/dt, y, dy/dt
    r_rk4 = r_init.copy()
    
    # Eq.8.33
    for t in tpoints_rk4:
        xpoints_rk4.append(r_rk4[0])
        vxpoints_rk4.append(r_rk4[1])
        ypoints_rk4.append(r_rk4[2])
        vypoints_rk4.append(r_rk4[3])
        k1 = h*rhs(r_rk4)  # all the k's are vectors
        k2 = h*rhs(r_rk4 + 0.5*k1)  # note: no explicit dependence on time of the RHSs
        k3 = h*rhs(r_rk4 + 0.5*k2)
        k4 = h*rhs(r_rk4 + k3)
        r_rk4 += (k1 + 2*k2 + 2*k3 + k4)/6

    return xpoints_rk4, vxpoints_rk4, ypoints_rk4, vypoints_rk4, tpoints_rk4


# The Adaptive method that would use at least 3 times RK4 method as described
# in the textbook pg.358-359
def adaptive_RK4(a,b,h,delta,r_init):
    """
    Parameters: a: Start time of simulation
                b: End time of simulation
                h: Initial time step
                delta: Target error-per-second
                r_init: Initial conditions in the order of [x, vx, y, vy]

    Returns: x_points_adaptive : Array that stores the x-position of particle
            vx_points_adaptive : Array that stores the x-velocity of particle
             y_points_adaptive : Array that stores the y-position of particle
            vy_points_adaptive : Array that stores the y-velocity of particle
              tpoints_adaptive : Array that stores time after each time step
                    h_adaptive : Araay that stores time step sizes
    """
    tpoints_adaptive = []
    x_points_adaptive = []
    vx_points_adaptive = []
    y_points_adaptive = []
    vy_points_adaptive = []
    h_adaptive = []
    
    t = a
    r_adaptive = r_init.copy()
    
    while t < b:
        x_points_adaptive.append(r_adaptive[0])
        vx_points_adaptive.append(r_adaptive[1])
        y_points_adaptive.append(r_adaptive[2])
        vy_points_adaptive.append(r_adaptive[3])
        
        tpoints_adaptive.append(t)
        h_adaptive.append(h)
        
        rho = 0
        
        while rho < 1:
            r_big_step = r_adaptive.copy()
            r_small_step = r_adaptive.copy()
            
            # Small step x(t+h) 2 times (x1)
            for _ in range(2):
                k1 = h * rhs(r_small_step)  # all the k's are vectors
                k2 = h * rhs(r_small_step + 0.5*k1)  # note: no explicit dependence on time of the RHSs
                k3 = h * rhs(r_small_step + 0.5*k2)
                k4 = h * rhs(r_small_step + k3)
                r_small_step += (k1 + 2*k2 + 2*k3 + k4)/6
            
            # Big step x(t+2h) (x2)
            k1 = 2*h * rhs(r_big_step)
            k2 = 2*h * rhs(r_big_step + 0.5*k1)
            k3 = 2*h * rhs(r_big_step + 0.5*k2)
            k4 = 2*h * rhs(r_big_step + k3)
            r_big_step += (k1 + 2*k2 + 2*k3 + k4)/6
            
            error_x = (r_small_step[0] - r_big_step[0])/30
            error_y = (r_small_step[2] - r_big_step[2])/30
            rho = h*delta / np.sqrt(error_x**2 + error_y**2)
            
            new_h = h * rho**0.25
            
            # Prevent new step size goes to infinity
            if new_h > 2*h:
                new_h = 2*h
            
            h = new_h
        
        # Calculate the next t location
        t = t + 2*h
        
        # Actual RK4 method using updated step (with 2h)
        for _ in range(2):
            k1 = h * rhs(r_adaptive)  # all the k's are vectors
            k2 = h * rhs(r_adaptive + 0.5*k1)  # note: no explicit dependence on time of the RHSs
            k3 = h * rhs(r_adaptive + 0.5*k2)
            k4 = h * rhs(r_adaptive + k3)
            r_adaptive += (k1 + 2*k2 + 2*k3 + k4)/6
            
    return x_points_adaptive, vx_points_adaptive, y_points_adaptive, vy_points_adaptive, \
            tpoints_adaptive, h_adaptive
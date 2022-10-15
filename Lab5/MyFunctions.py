# Import needed libraries
import numpy as np
def Euler_Cromer1D(t,x0,v0,m,k,dt):
    '''
    Performs the Euler cromer method to integrate
    the EOM of a relativistic particle in a spring. 
    Returns the position and velocity arrays of the particle
    '''
    from scipy.constants import c
    # Define zero vectors for mock variables
    u1 = np.zeros(len(t))
    u2 = np.zeros(len(t))
    u1dot = np.zeros(len(t))
    u2dot = np.zeros(len(t))

    # Set intial conditions
    u1[0] = x0
    u2[0] = v0

    # Run for loop from 1 to len(t)
    for i in range(1,len(t)):
        u2dot[i] = (-k/m)*u1[i-1]*(1-(u2[i-1]**2/c**2))**(3/2)
        u2[i] = u2[i-1] + dt*u2dot[i]

        u1dot[i] = u2[i]
        u1[i] = u1[i-1] + dt*u1dot[i]
    return u1
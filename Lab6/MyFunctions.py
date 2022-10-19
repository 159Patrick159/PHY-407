from re import A
import numpy as np

# Define acceleratino components from Q1a
def ax(x,y,a,b,sigma,epsilon,m):
    '''Returns the x-component of the acceleration
    for a particle subject to a Lennard-Jones
    potential. Where a and b are the x,y coordinates
    of the opposite particle'''
    term1 = (48*epsilon*sigma**12)/((x-a)**2 + (y-b)**2)**13
    term2 = (24*epsilon*sigma**6)/((x-a)**2 + (y-b)**2)**7
    return ((x-a)/m)*(term1 - term2)

def ay(x,y,a,b,sigma,epsilon,m):
    '''Returns the y-component of the acceleration
    for a particle subject to a Lennard-Jones
    potential.Where a and b are the x,y coordinates
    of the opposite particle'''
    term1 = (48*epsilon*sigma**12)/((x-a)**2 + (y-b)**2)**13
    term2 = (24*epsilon*sigma**6)/((x-a)**2 + (y-b)**2)**7
    return ((y-b)/m)*(term1 - term2)

def Verlet(t,dt,r01,r02):
    '''Performs the Verlet method to simulate
    the interaction of 2 particles under Lennard
    Jones potential. r01 and r02 are 1D arrays
    of length 2 where first entry is x0 and second
    entry is y0
    '''
    # Define constants for simulation
    epsilon = 1
    sigma = 1
    m = 1
    
    # Initialize array of zeros
    vx1 = np.zeros(len(t))
    vy1 = np.zeros(len(t))

    vx2 = np.zeros(len(t))
    vy2 = np.zeros(len(t))

    x1 = np.zeros(len(t))
    y1 = np.zeros(len(t))

    x2 = np.zeros(len(t))
    y2 = np.zeros(len(t))
    
    # Set initial conditions
    x1[0], y1[0] = r01[0], r01[1]
    x2[0], y2[0] = r02[0], r02[1]

    # Compute inital iteration for velocity
    vx1[1] = vx1[0] + dt/2*ax(r01[0],r01[1],r02[0],r02[1],sigma,epsilon,m)/2
    vy1[1] = vy1[0] + dt/2*ay(r01[0],r01[1],r02[0],r02[1],sigma,epsilon,m)/2

    vx2[1] = vx2[0] + dt/2*ax(r02[0],r02[1],r01[0],r01[1],sigma,epsilon,m)/2
    vy2[1] = vy2[0] + dt/2*ay(r02[0],r02[1],r01[0],r01[1],sigma,epsilon,m)/2  
    
    for i in range(len(t)-3):
        # Implement method for particle 1
        x1[i+2] = x1[i] + dt*vx1[i+1]
        y1[i+2] = y1[i] + dt*vy1[i+1]

        # Implement method for particle 2
        x2[i+2] = x2[i] + dt*vx2[i+1]
        y2[i+2] = y2[i] + dt*vy2[i+1]

        # Compute vector k = a(t+dt)
        kx1 = dt*ax(x1[i+2],y1[i+2],x2[i+2],y2[i+2],sigma,epsilon,m)/2
        ky1 = dt*ay(x1[i+2],y1[i+2],x2[i+2],y2[i+2],sigma,epsilon,m)/2

        # Compute vector k = a(t+dt)
        kx2 = dt*ax(x2[i+2],y2[i+2],x1[i+2],y1[i+2],sigma,epsilon,m)/2
        ky2 = dt*ay(x2[i+2],y2[i+2],x1[i+2],y1[i+2],sigma,epsilon,m)/2

        # Compute velocity v(t+dt)
        vx1[i+2] = vx1[i+1] + 0.5*kx1
        vy1[i+2] = vy1[i+1] + 0.5*ky1

        # Compute velocity v(t+dt)
        vx2[i+2] = vx2[i+1] + 0.5*kx2
        vy2[i+2] = vy2[i+1] + 0.5*ky2

        # Compute velocity v(t+3t/2)
        vx1[i+3] = vx1[i+1] + kx1
        vy1[i+3] = vy1[i+1] + ky1
    
        # Compute velocity v(t+3t/2)
        vx2[i+3] = vx2[i+1] + kx2
        vy2[i+3] = vy2[i+1] + ky2

    # Filter out the zeros in the array
    x1 = x1[np.arange(0,len(t),2)]
    y1 = y1[np.arange(0,len(t),2)]

    x2 = x2[np.arange(0,len(t),2)]
    y2 = y2[np.arange(0,len(t),2)]

    return(x1,y1,x2,y2)
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
    
    # Initialize arrays and set initial conditions
    vx1 = [0]
    vy1 = [0]

    vx2 = [0]
    vy2 = [0]

    x1 = [r01[0]]
    y1 = [r01[1]]

    x2 = [r02[0]]
    y2 = [r02[1]]


    # Compute first half-step velocity
    v0x1 = vx1[0] + dt/4*ax(x1[0],y1[0],x2[0],y2[0],sigma,epsilon,m) 
    v0y1 = vy1[0] + dt/4*ay(x1[0],y1[0],x2[0],y2[0],sigma,epsilon,m)

    v0x2 = vx2[0] + dt/4*ax(x2[0],y2[0],x1[0],y1[0],sigma,epsilon,m) 
    v0y2 = vy2[0] + dt/4*ay(x2[0],y2[0],x1[0],y1[0],sigma,epsilon,m)

    for i in range(len(t)-1):
        # Compute the next full step position for particle 1
        x1.append(x1[i] + dt*v0x1)
        y1.append(y1[i] + dt*v0y1)
        # Compute the next full step position for particle 2
        x2.append(x2[i] + dt*v0x2)
        y2.append(y2[i] + dt*v0y2)

        # Compute the k vector components for particle 1
        k1x = dt * ax(x1[i+1],y1[i+1],x2[i+1],y2[i+1],sigma,epsilon,m)/2
        k1y = dt * ay(x1[i+1],y1[i+1],x2[i+1],y2[i+1],sigma,epsilon,m)/2
        # Compute the k vector components for particle 2
        k2x = dt * ax(x2[i+1],y2[i+1],x1[i+1],y1[i+1],sigma,epsilon,m)/2
        k2y = dt * ay(x2[i+1],y2[i+1],x1[i+1],y1[i+1],sigma,epsilon,m)/2

        # Compute the next full step velocity components for particle 1
        vx1.append(v0x1 + 0.5*k1x)
        vy1.append(v0y1 + 0.5*k1y)
        # Compute the next full step velocity components for particle 2
        vx2.append(v0x2 + 0.5*k2x)
        vy2.append(v0y2 + 0.5*k2y)

        # Compute next half step velocity for particle 1
        v0x1 += k1x
        v0y1 += k1y
        # Compute next half step velocity for particle 2
        v0x2 += k2x
        v0y2 += k2y

    # Create arrays particle 1 and 2 with their positions and velocities
    r1 = np.array([x1,y1])
    r2 = np.array([x2,y2])
    v1 = np.array([vx1,vy1])
    v2 = np.array([vx2,vy2])
    return(r1,r2,v1,v2)

def Lennard_Jones(r):
    '''Returns the Lennar-Jones potential between two
    particles situated at (x1,y1) and (x2,y2), assuming 
    all the constants are equal to 1'''
    #r = np.sqrt((x1-x2)**2 + (y1-y2)**2)
    V = 4*(r**(-12)-r**(-6))
    return(V/2)

def KE(v):
    "Computes the kinetic energy of a particle of mass 1"
    return(0.5*v**2)

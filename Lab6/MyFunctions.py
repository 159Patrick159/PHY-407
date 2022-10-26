import numpy as np

############################### HEADER #######################################
# This file contains all the self-written function used for Lab6
##############################################################################

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

# Verlet method used for multi-body simulation
def Verlet_multibody(t,dt,r0,v0,N):
    """
    Parameters:  t: Time array
                dt: Time steps
                r0: [x0,y0], initial positions for all N particles
                v0: [vx0,vy0], initial velocities for all N particles
                 N: Number of particles simulating

    Returns: x_all_t: 2D numpy array containing x-pos for all N particles for all time t
             y_all_t: 2D numpy array containing y-pos for all N particles for all time t
           vx_all_t : 2D numpy array containing full-step x-vel for all N particles for all time t
           vy_all_t : 2D numpy array containing full-step y-vel for all N particles for all time t
       vxhalf_all_t : 2D numpy array containing half-step x-vel for all N particles for all time t
       vyhalf_all_t : 2D numpy array containing half-step y-vel for all N particles for all time t
                V_t : Potential energy of system as a function of time
                K_t : Kinetic energy of system as a function of time
    """
    # Define constants for simulation
    epsilon = 1
    sigma = 1
    m = 1
    
    # Full step position arrays for each particle as a function of t 
    #[Particle i][Time k]
    x_all_t = np.zeros((N, len(t)))
    y_all_t = np.zeros((N, len(t)))
    
    # Full step velocity arrays for each particle as a function of t
    #[Particle i][Time k]
    vx_all_t = np.zeros((N, len(t)))
    vy_all_t = np.zeros((N, len(t)))
    
    # Half-step velocity arrays for each particle as a function of t
    #[Particle i][Time k]
    vxhalf_all_t = np.zeros((N, len(t)))
    vyhalf_all_t = np.zeros((N, len(t)))
    
    # (Full-step) energy arrays for entire system as a function of t
    #[Time k]
    V_t = np.zeros(len(t))
    K_t = np.zeros(len(t))
    
    # Set initial conditions for all N particles from input
    x_all_t[:, 0] = r0[0]
    y_all_t[:, 0] = r0[1]
    vx_all_t[:, 0] = v0[0]
    vy_all_t[:, 0] = v0[1]
    
    # Calculate initial potential energy and kinetic energy
    for i in range(N):
        V_t[0] += potential_multibody(x_all_t[:,0], y_all_t[:,0], i)
        K_t[0] += KE(np.sqrt(vx_all_t[i][0]**2 + vy_all_t[i][0]**2))
    
    # Calculate the initial half-step velocity (Eq.8) and potential energy
    for i in range(N):
        ax_i, ay_i = acceleration_multibody(x_all_t[:,0], y_all_t[:,0], i, sigma,epsilon,m)
        vxhalf_all_t[i][0] = vx_all_t[i][0] + 0.5*dt*ax_i
        vyhalf_all_t[i][0] = vy_all_t[i][0] + 0.5*dt*ay_i
        
    
    # (Actual implementation of Verlet method) 
    # Looping over time step from 1 to T (shift index fwd once)
    for k in range(1,len(t)):
        
        # Calculate position of all particles (Eq.9)
        for i in range(N):
            x_all_t[i][k] = x_all_t[i][k-1] + dt * vxhalf_all_t[i][k-1]
            y_all_t[i][k] = y_all_t[i][k-1] + dt * vyhalf_all_t[i][k-1]
        
        # Update full-step and half-step velocities of all particles
        for i in range(N):
            ax_i_cur, ay_i_cur = acceleration_multibody(x_all_t[:,k], y_all_t[:,k], i, sigma,epsilon,m)
            
            # Compute the k vector components for particle i (Eq.10)
            kx_i = dt * ax_i_cur
            ky_i = dt * ay_i_cur
            
            # Compute the full-step velocity for particle i (Eq.11)
            vx_all_t[i][k] = vx_all_t[i][k-1] + 0.5 * kx_i
            vy_all_t[i][k] = vy_all_t[i][k-1] + 0.5 * ky_i
            
            # Compute the half-step velocity for particle i (Eq.12)
            vxhalf_all_t[i][k] = vxhalf_all_t[i][k-1] + kx_i
            vyhalf_all_t[i][k] = vyhalf_all_t[i][k-1] + ky_i
        
            # Calculate the total energy of the system at this time step
            V_t[k] += potential_multibody(x_all_t[:,k], y_all_t[:,k], i)
            K_t[k] += KE(np.sqrt(vx_all_t[i][k]**2 + vy_all_t[i][k]**2))
       
        #print(f"{100*(k+1)*dt/t[-1] : .2f}% done")
    return x_all_t, y_all_t, vx_all_t, vy_all_t, vxhalf_all_t, vyhalf_all_t, V_t, K_t
        
        
# Calculate Acceleration in multi-body simulation
def acceleration_multibody(x_all, y_all, i, sigma,epsilon,m):
    """
    Parameters: x_all : 1D numpy array containing x-pos for all N particles at set t-value
                y_all : 1D numpy array containing y-pos for all N particles at set t-value
                    i : The ith particle for the loop
    sigma, epsilon, m : Physical variables
    Returns: ax_i, ay_i: Acceleratioin x and y of the ith particles
    """
    ax_i = 0
    ay_i = 0
    for j in range(len(x_all)):
        if not (j==i):
            ax_i += ax(x_all[i], y_all[i], x_all[j], y_all[j], sigma,epsilon,m)
            ay_i += ay(x_all[i], y_all[i], x_all[j], y_all[j], sigma,epsilon,m)
    
    return (ax_i, ay_i)
            
# Calculate the potential energy per particle in multi-body simulation
def potential_multibody(x_all, y_all, i):
    """
    Parameters: x_all : 1D numpy array containing x-pos for all N particles at set t-value
                y_all : 1D numpy array containing y-pos for all N particles at set t-value
                    i : The ith particle for the loop
    Returns: total_potential_i : Potential energy of the ith particles
    """
    total_potential_i = 0
    for j in range(len(x_all)):
        if not (j==i):
            dx = x_all[j] - x_all[i]
            dy = y_all[j] - y_all[i]
            r = np.sqrt(dx**2 + dy**2)
            total_potential_i += Lennard_Jones(r)
            
    return total_potential_i
    

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

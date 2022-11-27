########################## 
# Contributors           # 
# Patrick Sandoval       #
# Kelvin Leong           #
##########################

################################# HEADER ##########################################
# This python script holds the code for question 1 of lab 10 for PHY407           #  
# all plots of random walk are also generated by this script. Some auxiliary      #
# functions are imported from the script MyFunctions.py                           #
###################################################################################

####################################### Q1a #######################################
# Import needed libraries
import numpy as np
import matplotlib.pyplot as plt
from pylab import clf, pause, draw

# Import direction function
from MyFunctions import nextmove

# Define domain with dimensions of box
L = 101
Z = np.zeros((L+1,L+1))
# Set inital position
x0, y0 = Z.shape[0]//2, Z.shape[1]//2
Z[y0,x0] = 1

# Store coordinates for plots
xstore = [x0]
ystore = [y0]
# Generate array for plotting animation
# xdomain = np.arange(0,L+1,1)
# ydomain = np.arange(0,L+1,1)
# Plot initial state
# plt.xlim([0,L])
# plt.ylim([0,L])
# plt.contourf(xdomain,ydomain,Z,cmap='inferno_r',origin="lower")
for i in range(500):
    # Generate direcion
    x,y = nextmove(x0,y0)
    # Check if new direction is within bounds
    if np.logical_and(np.logical_and(x>=0,x<=L),np.logical_and(y>=0,y<=L)):
        # Apply new diection
        Z[y,x] = 1
        # Set old position to zeor
        Z[y0,x0] = 0
        xstore.append(x)
        ystore.append(y)
        # Animation uncomment if curious
        # clf()
        # plt.contourf(xdomain,ydomain,Z,cmap='inferno',origin="lower")
        # plt.xlim([0,L])
        # plt.ylim([0,L])
        # plt.xlabel("X (Spatial)")
        # plt.ylabel("Y (Spatial)")
        # plt.title(f"Motion of Particle in {L}x{L} Lattice")
        # draw()
        # pause(0.001)
        # Set new coords to old coords
        x0, y0 = x, y

# Generate time array
t = np.linspace(0,500,len(xstore))
fig, (a0,a1,a2) = plt.subplots(figsize=(12,4),ncols=3)
a0.plot(t,xstore,c='k')
a0.grid(ls='--')
a0.set_title("X-Component of Motion",fontsize=16)
a0.set_xlabel("Time",fontsize=14)
a0.set_ylabel("X Displacement",fontsize=14)

a1.plot(t,ystore,c='r')
a1.grid(ls='--')
a1.set_title("Y-Component of Motion",fontsize=14)
a1.set_xlabel("Time",fontsize=14)
a1.set_ylabel("Y Displacement",fontsize=14)

a2.plot(xstore,ystore,c='g')
a2.grid(ls='--')
a2.set_title("2-Dimensional Motion",fontsize=16)
a2.set_ylabel("Y Displacement",fontsize=14)
a2.set_xlabel("X Displacement",fontsize=14)

plt.tight_layout()
#plt.savefig("Lab10\Q1aPlot.pdf")
plt.show()


####################################### Q1b #######################################
import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import nextmove, Collision

# Generate grid with edges of ones
L = 101
grid = np.ones((L+3,L+3))
grid[1:grid.shape[1]-1,1:grid.shape[0]-1] = 0

# Set initial particle in center of grid
x0, y0 = grid.shape[0]//2, grid.shape[1]//2
grid[y0,x0] = 1
end = False
while not end:
    # Check if there was collision at previous step
    while not Collision(x0,y0,grid):
        # If are no collisions at x0 and y0 we continue RW
        x,y = nextmove(x0,y0)
        #Update particle position
        grid[y,x] = 1
        # Remove previous instance
        grid[y0,x0] = 0
        x0,y0 = x,y
    # Check for collision at origin
    if Collision( grid.shape[0]//2, grid.shape[1]//2,grid):
        end = True
    # If there is a collision we reset particle to center
    x0, y0 = grid.shape[0]//2, grid.shape[1]//2

xdom = np.arange(0,L+1,1)
ydom = np.arange(0,L+1,1)
grid_plot = grid[1:grid.shape[1]-1,1:grid.shape[0]-1]
plt.contourf(xdom,ydom,grid_plot,cmap="inferno",origin='lower')
plt.xlabel("X (Spatial)",fontsize=14)
plt.ylabel("Y (Spatial)",fontsize=14)
plt.title("Diffusion-Limited Aggregation Model",fontsize=16)
plt.grid(ls='--')
plt.tight_layout()
#plt.savefig("Lab10/Q1bPlot.pdf")
plt.show()
        

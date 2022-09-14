# Collaborators:
# Patrick Sandoval
# Kelvin Leong

########################### HEADER ###############################
# This python script times the matrix multiplication of 2 NxN matrices
# The code will print out the time it takes to compute the operation and
# it will plot the time as a function of N

############################################ Q3 ################################################
# Import needed libraries
import numpy as np
import matplotlib.pyplot as plt
from time import time
times = []
Ntop = 120
# Define A and B based on N
for N in range(2,Ntop):
    A = np.ones([N,N],float)*3
    B = np.ones([N,N],float)

    # Use the code snippet to compute new matric C
    C = np.zeros([N,N],float)
    # Star timing 
    start = time()
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i,j] += A[i,k]*B[k,j]
    # End timing and compute time ellapsed
    end = time()
    diff = end-start
    times.append(diff)

# Use numpy.dot to perform similar operation
np_times = []
for N in range(2,Ntop):
    A = np.ones([N,N],float)*3
    B = np.ones([N,N],float)
    start = time()
    C = np.dot(A,B)
    end = time()
    diff = end-start
    np_times.append(diff)

# Create domain for plotting
Ndomain = np.arange(2,Ntop,1)

fig, ((a0,a1),(a2,a3)) = plt.subplots(figsize=(10,5),ncols=2,nrows=2,gridspec_kw={'height_ratios': [3, 1]})
a0.scatter(Ndomain,times,c='k',s=6)
a0.grid(ls='--')
a0.set_title("Timing matrix multiplication for N",fontsize=16)
#a0.set_xlabel("Size N",fontsize=14)
a0.set_ylabel("Time [s]",fontsize=14)

a1.scatter(Ndomain**3,times,c='r',s=6)
a1.grid(ls='--')
a1.set_title(r"Timing matrix multiplication for N$^3$",fontsize=16)
#a1.set_xlabel(r"N$^3$",fontsize=14)
a1.set_ylabel("Time [s]",fontsize=14)

a2.scatter(Ndomain,np_times,c='gray',s=6)
a2.grid(ls='--')
a2.set_title("Numpy.dot Performance",fontsize=16)
a2.set_xlabel("Size N",fontsize=14)
a2.set_ylabel("Time [s]",fontsize=14)

a3.scatter(Ndomain**3,np_times,c='orange',s=6)
a3.grid(ls='--')
a3.set_title(r"Numpy.dot Performance",fontsize=16)
a3.set_xlabel(r"N$^3$",fontsize=14)
a3.set_ylabel("Time [s]",fontsize=14)

plt.tight_layout()
plt.savefig("Q3-Plot.png")
plt.show()




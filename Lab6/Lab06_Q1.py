# Contributors: 
# Patrick Sadnoval
# Kelvin Leong


####################################### HEADER ###########################################
# This python script holds the code, pseudo-code for the simulations for the 
# 2-dimensional molecular simulation for Q1 and all its subsections
######################################## Q1.a ############################################
import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import Verlet

# Define time step 
dt = 0.01
t = np.arange(0,100*dt,dt)

r11 = [4,4]
r21 = [5.2,4]

r12 = [4.5,4]
r22 = [5.2,4]

r13 = [2,3]
r23 = [3.5,4.4]

r1,r2 = Verlet(t,dt,r11,r21)
r3,r4 = Verlet(t,dt,r12,r22)
r5,r6 = Verlet(t,dt,r13,r23)

fig, (a0,a1,a2) = plt.subplots(figsize=(12,4),ncols=3)
a0.plot(r1[0],r1[1],'.',label="Particle 1",c='k')
a0.plot(r2[0],r2[1],'.',label="Particle 2",c='r')
a0.legend()

a1.plot(r3[0],r3[1],'.',label="Particle 1",c='k')
a1.plot(r4[0],r4[1],'.',label="Particle 2",c='r')
a1.legend()

a2.plot(r5[0],r5[1],'.',label="Particle 1",c='k')
a2.plot(r6[0],r6[1],'.',label="Particle 2",c='r')
a2.legend()

plt.tight_layout()
plt.savefig("Q1a.pdf")
plt.show()
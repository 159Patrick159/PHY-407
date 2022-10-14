# Contritbutors
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code, pseudo-code and plots
# for Q3 and all its subsections for Lab5
######################################## Q3.a ############################################
# Import needed libraries
from sqlite3 import Time
import numpy as np
import matplotlib.pyplot as plt

# Load data fiels
SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')

# Plot original data
plt.contourf(Longitude, Times, SLP)
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa)')
plt.colorbar()
plt.tight_layout()
plt.savefig("Q3a.pdf")
plt.show()


# Take 2D fourier transform and center zeroth frequency
G2 = np.fft.fftshift(np.fft.fft2(SLP))

# Compute the step of Longitude and time
dL = (Longitude[-1] - Longitude[0])/len(Longitude)
dt = (Times[-1] - Times[0])/len(Times)

u = np.fft.fftshift(np.fft.fftfreq(len(Longitude),dL))
v = np.fft.fftshift(np.fft.fftfreq(len(Times),dt))

plt.contourf(u,v,np.abs(G2))
#plt.contourf(np.abs(G2))
plt.xlabel(r"1/$\theta$")
plt.ylabel(r"1/s")
plt.colorbar()
plt.tight_layout()
plt.savefig("Q3aFourier.pdf")
plt.show()

# Set all other wave numbers fromm = 3 and m = 5 to zero
G23 = np.copy(G2)
G25 = np.copy(G2)
m = 3
n = 5
G23[:,:len(u)//2 + m] = 0
G23[:,len(u)//2 + m+1:] = 0

G25[:,:len(u)//2 + n] = 0
G25[:,len(u)//2 + n+1:] = 0

# Take the inverse fourier of cleaned spectrum
f23 = np.fft.ifft2(G23)
f25 = np.fft.ifft2(G25)

plt.contourf(Longitude,Times,np.abs(f23))
plt.colorbar()
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa) for m=3')
plt.savefig("Q3a3.pdf")
plt.show()

plt.contourf(Longitude,Times,np.abs(f25))
plt.colorbar()
plt.xlabel('longitude(degrees)')
plt.ylabel('days since Jan. 1 2015')
plt.title('SLP anomaly (hPa) for m=5')
plt.savefig("Q3a5.pdf")
plt.show()



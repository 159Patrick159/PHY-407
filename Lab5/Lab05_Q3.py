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
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Load data fiels
SLP = np.loadtxt('SLP.txt')
Longitude = np.loadtxt('lon.txt')
Times = np.loadtxt('times.txt')

# Plot original data
fig, a0 = plt.subplots(figsize=(8,6))
im0 = a0.contourf(Longitude, Times, SLP)
a0.set_xlabel(r'Longitude($^{\circ}$)',fontsize=18)
a0.set_ylabel('Days since Jan. 1 2015',fontsize=18)
a0.set_title('SLP anomaly (hPa)',fontsize=20)
divider = make_axes_locatable(a0)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im0, cax=cax, orientation='vertical',label="Pressure (hPa)")
plt.tight_layout()
plt.savefig("Q3a.pdf")


# Take 2D fourier transform and center zeroth frequency
G2 = np.fft.fft2(SLP)

# Compute the step of Longitude and time
dL = (Longitude[-1] - Longitude[0])/len(Longitude)
dt = (Times[-1] - Times[0])/len(Times)

u = np.fft.fftshift(np.fft.fftfreq(len(Longitude),dL))
v = np.fft.fftshift(np.fft.fftfreq(len(Times),dt))

fig, a0 = plt.subplots(figsize=(8,6))
im0 = a0.contourf(u,v,np.abs(np.fft.fftshift(G2)))
a0.set_xlabel(r"Cycles/$\theta$",fontsize=18)
a0.set_ylabel(r"Cycles/Days",fontsize=18)
a0.set_title(r"|$\bar{G}$(u,v)|$^2$",fontsize=20)
divider = make_axes_locatable(a0)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im0, cax=cax, orientation='vertical')
plt.tight_layout()
plt.savefig("Q3aFourier.pdf")

# Set all other wave numbers fromm = 3 and m = 5 to zero
G23 = np.copy(G2)
G25 = np.copy(G2)

m = 3
n = 5
G23[:,:m] = 0
G23[:,m+1:] = 0

G25[:,:n] = 0
G25[:,n+1:] = 0

# Take the inverse fourier of cleaned spectrum
f23 = np.fft.ifft2(G23)
f25 = np.fft.ifft2(G25)

fig, a0 = plt.subplots(figsize=(8,6))
im0 = a0.contourf(Longitude,Times,f23.real)
a0.set_xlabel('Longitude(degrees)',fontsize=18)
a0.set_ylabel('Days since Jan. 1 2015',fontsize=18)
a0.set_title(f'SLP anomaly (hPa) for m={m}',fontsize=20)
divider = make_axes_locatable(a0)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im0, cax=cax, orientation='vertical')
plt.tight_layout()
plt.savefig("Q3a3.pdf")

fig, a0 = plt.subplots(figsize=(8,6))
im0 = a0.contourf(Longitude,Times,f25.real)
a0.set_xlabel('Longitude(degrees)',fontsize=18)
a0.set_ylabel('Days since Jan. 1 2015',fontsize=18)
a0.set_title(f'SLP anomaly (hPa) for m={m}',fontsize=20)
divider = make_axes_locatable(a0)
cax = divider.append_axes('right', size='5%', pad=0.05)
fig.colorbar(im0, cax=cax, orientation='vertical')
plt.tight_layout()
plt.savefig("Q3a5.pdf")




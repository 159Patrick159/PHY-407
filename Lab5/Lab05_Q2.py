# Contritbutors
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code, pseudo-code and plots
# for Q2 and all its subsections for Lab5
##############################################################################
# Pseudocode: 
# Import the necessary libraries and functions
# (Q2.b,c)
# Read data from GraviteaTime.wav into 2 channels (in units of Hz)
# Use the number of samples in either channel to create a time domain
# Plot the data of both channels on the time domain, and limit/zoom-in the 
# graph in the range of 0.02s to 0.05s
# (Q2.d)
# Perform Fourier transform on both the channels to the frequency domain and 
# plot their amplitude
# Set the transformed data who is greater than 880 Hz to be zero, and plot
# their (filtered) amplitude
# Perform inverse Fourier transform on both filter data back to the time domain,
# and plot their filtered time series (in the same range of 0.02s to 0.05s)
# (Q2.e)
# Write both the filtered time series channels to a new file called 
# GraviteaTime lpf.wav
######################################## Q2.b ############################################
# Import needed functions
from scipy.io.wavfile import read, write
from numpy import empty
import numpy as np
import matplotlib.pyplot as plt

# read the data into two stereo channels
# sample is the sampling rate, data is the data in each channel
# dimensions [2, nsamples]
sample, data = read('GraviteaTime.wav')
# sample is the sampling frequency, 44100 Hz
# separate into channels
channel_0 = data[:, 0]
channel_1 = data[:, 1]
N_Points = channel_0.size

# Create time domain
dt = 1/sample
T = len(channel_0)*dt
t = np.arange(0,T,dt)

fig, (a0,a1) = plt.subplots(figsize=(9,4),nrows=2)

a0.plot(t,channel_0,c='k')
a0.set_title("Channel 0")
a0.set_xlabel("Time [s]")
a0.set_xlim([0.02,0.05])

a1.plot(t,channel_1,c='k')
a1.set_title("Channel 1")
a1.set_xlabel("Time [s]")
a1.set_xlim([0.02,0.05])
plt.tight_layout()
plt.savefig("Q2a.pdf")
#%%
######################################## Q2.d ############################################
# Take fourier transform of each individual channel
# And shift to zero frequency
G0 = np.fft.fftshift(np.fft.fft(channel_0))
G1 = np.fft.fftshift(np.fft.fft(channel_1))

# Create frequency domain
f = np.fft.fftshift(np.fft.fftfreq(len(t),dt))

# Plot original fourier spectrum
fig, (a0,a1) = plt.subplots(figsize=(9,4),nrows=2)
a0.plot(f,np.abs(G0),c='k')
a0.set_xlabel("Frequency [Hz]")
a0.set_ylabel(r"|$\bar{G}(f)$|$^2$")
a0.set_title("Channel 0 Fourier Spectrum")
a0.set_xlim([np.min(f)/2,np.max(f)/2])

a1.plot(f,np.abs(G1),c='k')
a1.set_xlabel("Frequency [Hz]")
a1.set_ylabel(r"|$\bar{G}(f)$|$^2$")
a1.set_title("Channel 1 Fourier Spectrum")
a1.set_xlim([np.min(f)/2,np.max(f)/2])

plt.tight_layout()
plt.savefig("Q2Fourier_Spec.pdf")

#%%
# Set values to zero for Freq>880Hz 
G0[np.where(abs(f)>880)] = 0
G1[np.where(abs(f)>880)] = 0

# Plot modified fourier spectrum
fig, (a0,a1) = plt.subplots(figsize=(9,4),nrows=2)
a0.plot(f,np.abs(G0),c='k')
a0.set_xlabel("Frequency [Hz]")
a0.set_ylabel(r"|$\bar{G}(f)$|$^2$")
a0.set_title("Channel 0 Modified Fourier Spectrum")
a0.set_xlim([np.min(f)/2,np.max(f)/2])

a1.plot(f,np.abs(G1),c='k')
a1.set_xlabel("Frequency [Hz]")
a1.set_ylabel(r"|$\bar{G}(f)$|$^2$")
a1.set_title("Channel 1 Modified Fourier Spectrum")
a1.set_xlim([np.min(f)/2,np.max(f)/2])

plt.tight_layout()
plt.savefig("Q2ModFourier_Spec.pdf")
#%%
# Take inverse fourier transfer of new Gs
mod_channel0 = np.fft.ifft(np.fft.ifftshift(G0))
mod_channel1 = np.fft.ifft(np.fft.ifftshift(G1))

fig, (a0,a1) = plt.subplots(figsize=(9,4),nrows=2)

a0.plot(t,mod_channel0.real,c='k')
a0.set_title("Channel 0")
a0.set_xlabel("Time [s]")
a0.set_xlim([0.02,0.05])

a1.plot(t,mod_channel1.real,c='k')
a1.set_title("Channel 1")
a1.set_xlabel("Time [s]")
a1.set_xlim([0.02,0.05])
plt.tight_layout()
plt.savefig("Q2d.pdf")

#%%
################################## Q2e #######################################
# this creates an empty array data_out with the same shape as "data"
# (2 x N_Points) and the same type as "data" (int16)
# data_out = empty(data.shape, dtype = data.dtype)
# # fill data_out
data_out = empty(data.shape, dtype = data.dtype)
data_out[:, 0] = mod_channel0.real
data_out[:, 1] = mod_channel1.real
write('GraviteaTime_lpf.wav', sample, data_out)
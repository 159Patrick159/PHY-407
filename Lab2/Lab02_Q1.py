# Collaborators:
# Patrick Sandoval
# Kelvin Leong

####################################### HEADER ###########################################
# This python script holds the code and pseudocode for all subquestion in Q1 of Lab02.
# Additionally we explain some parts of the code and provide useful comments. The dataset
# must be located in the same working directory of this script for the code to run properly.


######################################## Q1.a ############################################
# Dynamically collect the current working directory with os.getcwd()
# Add the dataset filename to the cwd string
# Call np.loadtxt() on the created filename to load the data
#
# Equation (1) does not have the possibility of taking the square root of a negative value
# therefore we do not need a checker when implementing this algorithm
# 
# We define a function for computing the mean of a dataset of size N
# Import function into this script
# Call the function and store the mean of the datase in a variable
# Define variable sum
# Define variable for data length
# Run a for loop for N iterations
#   Compute the square of the difference between the ith element and the mean
#   Add the result of operation to sum variable 
# Complete the rest of the fomula (1)
#
# Equation (2) does have the possibilit of taking the square root of a negative value
# so we will need to implement a checker to stop the code if this indeed happens
# 
# Define variable for square sum
# Define variable for normal sum
# Run a for loop for N iterations
#   square the ith entry and add value to square sum variable
#   add ith value like a normal sum
# Compute square mean by dividing normal sum by N and squaring the result
# Compute term square sum - N*square mean
# Check if difference is >0
# If condition fails print error message
# If condition passes then we multiply by 1/(N-1) and take the square root
# 
# Compute the "true" mean with numpy function
# Compute and print the relative errors for each method

######################################## Q1.b ############################################
# Import needed libraries
from logging import raiseExceptions
import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

# Dynamically get current working directory
current_directory = os.getcwd()
# Add filename to directory path
file_path = current_directory + "\cdata.txt"
# Load data
data = np.loadtxt(file_path)
N = len(data)

# Method 1
from MyFunctions import find_mean
xbar = find_mean(data)
sum1 = 0
for val in data:
    sum1 += (val - xbar)**2
std1 = np.sqrt((1/(N-1))*sum1)

# Method 2
square_sum = 0
norm_sum = 0
for val in data:
    square_sum += val**2
    norm_sum += val

square_xbar = (norm_sum/N)**2
term = square_sum-(N*square_xbar)

# Implement checker
flag=True
if term<0:
    flag = False
else:
    std2 = np.sqrt((1/(N-1))*term)

#Compute relative error
true = np.std(data,ddof=1)
rel_err1 = (std1-true)/true
if flag:
    rel_err2 = (std2-true)/true

# Print relative errors
print("Two Pass Relative Error:",rel_err1)
if flag:
    print("One Pass Relative Error:",rel_err2)
else:
    print("One Pass: Cannot take sqrt of negative numbers") 

print()

######################################## Q1.c ############################################

# Define mean, sigma and lenght for both Gaussians
mean1, sigma1, N1 = 0.,1.,2000
mean2, sigma2, N2 = 1e7, 1., 2000

# Generate both Gaussian distributions for given parameters
dist1 = np.random.normal(mean1,sigma1,N1)
dist2 = np.random.normal(mean2,sigma2,N2)

# Define tick direction 
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'

fig, (a0,a1) = plt.subplots(figsize=(10,4),ncols=2)
a0.hist(dist1,bins = 'fd',density=True, color = 'k', histtype = 'step', alpha = 0.5)
a0.axvline(x=np.mean(dist1),ls="--",c='r',label=r"$\mu$")
a0.axvline(x=np.mean(dist1)+np.std(dist1),ls="--",c='k',label=r"$\sigma$")
a0.axvline(x=np.mean(dist1)-np.std(dist1),ls="--",c='k')
a0.plot(dist1,np.full_like(dist1, 0.01), '|k', markeredgewidth = 1, alpha = 0.1)
a0.xaxis.set_minor_locator(MultipleLocator(0.2))
a0.yaxis.set_minor_locator(MultipleLocator(0.01))
a0.legend()
a0.set_title(r"Gaussian Distribution with $\mu$=0")

a1.hist(dist2,bins = 'fd',density=True, color = 'k', histtype = 'step', alpha = 0.5)
a1.axvline(x=np.mean(dist2),ls="--",c='r',label=r"$\mu$")
a1.axvline(x=np.mean(dist2)+np.std(dist2),ls="--",c='k',label=r"$\sigma$")
a1.axvline(x=np.mean(dist2)-np.std(dist2),ls="--",c='k')
a1.plot(dist2,np.full_like(dist1, 0.01), '|k', markeredgewidth = 1, alpha = 0.1)
a1.xaxis.set_minor_locator(MultipleLocator(0.2))
a1.yaxis.set_minor_locator(MultipleLocator(0.01))
a1.legend()
a1.set_title(r"Gaussian Distribution with $\mu=10^{7}$")
plt.tight_layout()
plt.savefig("Q1cPlot.png")
plt.show()

# Compute the standard deviations using eq.1 an eq.2 from lab manual
# For dist 1
xbar11 = find_mean(dist1)
sum11 = 0
for val in dist1:
    sum11 += (val - xbar11)**2
std11 = np.sqrt((1/(N1-1))*sum11)

square_sum21 = 0
norm_sum21 = 0
for val in dist1:
    square_sum21 += val**2
    norm_sum21 += val

square_xbar21 = (norm_sum21/N1)**2
term21 = square_sum21-(N1*square_xbar21)

# Implement checker
flag1 = True
if term21<0:
    flag1 = False
else:
    std21 = np.sqrt((1/(N1-1))*term21)

# For dist 2
xbar12 = find_mean(dist2)
sum12 = 0
for val in data:
    sum12 += (val - xbar12)**2
std12 = np.sqrt((1/(N2-1))*sum12)

square_sum22 = 0
norm_sum22 = 0
for val in dist2:
    square_sum22 += val**2
    norm_sum22 += val

square_xbar22 = (norm_sum22/N2)**2
term22 = square_sum22-(N2*square_xbar22)

# Implement checker
flag2 = True
if term22<0:
    flag2 = False
else:
    std22 = np.sqrt((1/(N2-1))*term22)


# Compute relative error
# For dist 1
rel11 = (std11 - sigma1)/(sigma1)
if flag1:
    rel12 = (std21 - sigma2)/sigma2
else:
    print("Negative sqrt in dist1 method 2")
    raiseExceptions("No imaginary numbers")

# For dis 2
rel21 = (std21 - sigma2)/(sigma2)
if flag2:
    rel22 = (std22 - sigma2)/(sigma2)
else:
    print("Negative sqrt in dist2 method 2")
# Print relative errors
print("Distribution with mean=0:")
print("Relative error for two pass:",rel11)
if flag1:
    print("Relative error for one pass:",rel12)
print()
print("Distribution with mean=1e7:")
print("Relative error for two pass:",rel21)
if flag2:
    print("Relative error for one pass:",rel22)
print()
######################################## Q1.d ############################################
# 
# Call find_mean to compute the mean of the data set
# Square the mean and multiply by N
# Compute the difference term using np.sum() and specifying the input to be the square of 
# the data entries.
# Create condition for negative values
# Compute relative errors
# Print out result

# Compute mean
#xbar3 = find_mean(data)
# Compute term of interest in a single line
term3 = np.sum(data**2) - N*(find_mean(data))**2

# Implement checker
flag=True
if term3<0:
    flag = False
else:
    std2 = np.sqrt((1/(N-1))*term3)

#Compute relative error
true = np.std(data,ddof=1)

if flag:
    rel_err2 = (std2-true)/true

# Print relative errors
if flag:
    print("One Pass Improved Relative Error:",rel_err2)
else:
    print("One Pass: Cannot take sqrt of negative numbers") 

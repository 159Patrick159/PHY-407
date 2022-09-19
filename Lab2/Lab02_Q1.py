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
# Define variable for sum
# Run a for loop for N iterations
#   square the ith entry and subtract the squared mean times N from it
#   add value of operation to sum
# Once for loop is done check if sum variable is >0
# If condition fails print error message
# If condition passes then we multiply by 1/(N-1) and take the square root
# 
# Compute the "true" mean with numpy function
# Compute and print the relative errors for each method

######################################## Q1.b ############################################
# Import needed libraries
import numpy as np
import matplotlib.pyplot as plt
import os

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
    sum1 =+ (val - xbar)**2
std1 = np.sqrt((1/(N-1))*sum1)

# Method 2
#sum2 = 0
#for val in data:
    #sum2 =+ val**2 - N*xbar**2
# Create checker for negative sum

# flag = True
# if sum2 < 0:
#     flag = False
# else:
#     std2 = np.sqrt((1/(N-1))*sum2)

# New approach
term = np.sum(data**2) - N*xbar**2

# Impement checker
flag=True
if term<0:
    flag = False
else:
    std2 = np.sqrt((1/(N-1))*term)

#Compute reative error
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



######################################## Q1.c ############################################

# Define mean, sigma and lenght for both Gaussians
mean1, sigma1, N1 = 0.,1.,2000
mean2, sigma2, N2 = 1e7, 1., 2000

# Generate both Gaussian distributions for given parameters
dist1 = np.random.normal(mean1,sigma1,N1)
dist2 = np.random.normal(mean2,sigma2,N2)

# Compute the standard deviations using np.std
std1 = np.std(dist1,ddof=1)
std2 = np.std(dist2,ddof=1)

# Compute relative error
rel1 = (std1-sigma1)/sigma1
rel2 = (std2-sigma2)/sigma2

# Print relative errors
print("Relative error for mean = 0:",rel1)
print("Relative error for mean = 1e7:",rel2)
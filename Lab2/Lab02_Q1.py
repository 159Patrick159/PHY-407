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
# Run a for loop for N iterations
#   square the ith entry and add value to square sum variable
# Compute square mean through square sum
# Compute term square sum - N*square mean
# Check if difference is >0
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
square_sum = 0
for val in data:
    square_sum =+ val**2
square_xbar = square_sum/N**2
term = square_sum-N*square_xbar
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

print()

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
print()
######################################## Q1.d ############################################
# The issue with eq.2 from the lab manual is that sum of the squared elements loses some
# of its decimal information due to the computers roundoff error of e16 sif figs. Because
# this quantity is really big the roundoff error make the quantity
# lose important decimal value information. And same thing happens with the term
# N*xbar**2 where roundouff error reduces the precision of the value, and when we take the
# difference of these two quantites the precision of the difference is no longer e16. 
#
# So we need to eliminate the 1 pass element and compute the difference term in a single run
# 
# Call find_mean to compute the mean of the data set
# Compute the difference term using np.sum() and specifying the input to be the square of 
# the data entries.
# Create condition for negative values
# Compute relative errors
# Print out result

# Compute mean
xbar3 = find_mean(data)
# Compute term of interest in a single line
term3 = np.sum(data**2) - N*xbar3**2

# Impement checker
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

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 10:56:23 2022

@Collaborator:  Kelvin Leong
                Patrick Sandoval
"""
############################ HEADER ##########################################
# This python code contains the explanation, codes and plots for calculating
# numerical roundoff errors in the Q4 exercise
##############################################################################
# Import the necessary libraries
import numpy as np
import matplotlib.pyplot as plt
import MyFunctions as myf

#%%
################################# Q4a #######################################
# Pick 500 points in the range 0.98 < u < 1.02
N = 500
u_range = np.linspace(0.98, 1.02, N)

# Calculate p(u) and q(u)
p_vals = myf.p(u_range)
q_vals = myf.q(u_range)

# Plot p(u) and q(u) vs. u
plt.figure(figsize=(8,8))
plt.plot(u_range, p_vals, 'r.', label='p(u)')
plt.plot(u_range, q_vals, 'k.', label='q(u)')
plt.xlabel("u value", fontsize=16)
plt.ylabel("Polynomial numerical value", fontsize=16)
plt.title("Comparison of numerical value between expanded and simplified \nform of polynomial", fontsize=14)
plt.grid('on')
plt.legend()
plt.savefig("Q4a.png")

#%%
################################# Q4b #######################################
# Calculate the difference p(u)-q(u)
diff_vals = p_vals - q_vals

# Plot the difference p(u)-q(u)
plt.figure(figsize=(8,8))
plt.plot(u_range, diff_vals, ".", label='p(u)-q(u)')
plt.xlabel("u value", fontsize=16)
plt.ylabel("p(u)-q(u) numerical value", fontsize=16)
plt.title("Difference in numerical value between expanded and simplified \nform of polynomial", fontsize=14)
plt.grid('on')
plt.legend()
plt.savefig("Q4b-Plot.png")

# Plot the histogram of the difference p(u)-q(u) 
bin_array = np.linspace(min(diff_vals), max(diff_vals), 20)
plt.figure(figsize=(8,8), tight_layout=True)
plt.hist(diff_vals, bin_array, align='left', label='Data')
plt.xlabel("p(u)-q(u)",fontsize=16)
plt.ylabel("Number of occurance",fontsize=16)
plt.title("Distribution of numerical value of p(u)-q(u)",fontsize=16)
plt.savefig("Q4b-Histogram.png")

# Calculate the standard deviation of the histogram (distribution)
std_dev = np.std(diff_vals, ddof=1)


# Calculate the roundoff error using equation 3 in the lab (C=1e-16)
rdf_error = myf.q4b_roundoff_error(u_range)
# Take the median of the roundoff errors over the u values near 1
mean_rdf_error = np.mean(rdf_error)

# Comparison between staandard deviation of distribution and the roundoff error estimate
print(f"Standard deviation of distribution: {std_dev: .3e}")
print(f"Roundoff error estimate: {mean_rdf_error: .3e}")

#%%
############################### Q4c #########################################
# Pick points greater than 0.98 but less than 1.0, calculate p(u) and q(u)
u_range = np.linspace(0.980, 0.985, 500)
p_vals = myf.p(u_range)
q_vals = myf.q(u_range)

# Calculate and Plot this absolute error
err = abs(p_vals-q_vals)/abs(p_vals)

plt.figure(figsize=(8,8))
plt.plot(u_range, err, ".")
plt.xlabel("u value",fontsize=16)
plt.ylabel("|p-q|/|p|",fontsize=16)
plt.title("|p-q|/|p| for 0.980 < u < 0.985", fontsize=16)
#plt.yscale('log')
plt.grid('on')
plt.savefig("Q4c-abs.png")

#%%
u_range = np.arange(0.982, 0.984, 0.0001) #0.9831

# Calculate fractional error of p(u)-q(u) and take the absolute
fractional_error = myf.q4c_fractional_error(u_range)
abs_fractional_error = abs(fractional_error)

# Plot the absolute fractional error
plt.figure(figsize=(8,8))
plt.plot(u_range, abs_fractional_error*100, ".", markersize=10)
plt.xlabel("u value",fontsize=16)
plt.ylabel("Absolute fractional error (%)",fontsize=16)
plt.title("Absolute fractional error for 0.982 < u < 0.984", fontsize=16)
#plt.yscale('log')
plt.ylim([0,500])
plt.grid('on')
plt.savefig("Q4c-frac-error.png")

print(abs_fractional_error)
print(f"Median of the absolute fractional error: {np.median(abs_fractional_error)*100:.2f} %")


#%%
############################## Q4d ##########################################

# Pick 500 points in the range 0.98 < u < 1.02, calculate f(u), f(u)-1
u_range = np.linspace(0.98, 1.02, 500)
f_vals = myf.q4d_f(u_range)
product_quotient_diff = f_vals - 1

# Plotting numerical value of f(u)-1
plt.figure(figsize=(8,8))
plt.plot(u_range, product_quotient_diff, ".", label='f-1')
plt.xlabel("u value", fontsize=16)
plt.ylabel("f(u) - 1 = u^8/(u^4 * u^4) - 1", fontsize=16)
plt.title("Numerical value for f(u)-1", fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig("Q4d.png")

# Calculate the standard deviation of f(u)-1
std_dev = np.std(product_quotient_diff, ddof=1)

# Calculate roundoff error according to eq4.5 (ignoring sqrt(2))
quotient_rdoff_error = u_range**8 / (u_range**4 * u_range**4)
C = 1e-16
roundOffError = np.mean(C * quotient_rdoff_error )

print(f"Standard deviation of f(u) distribution: {std_dev}")
print(f"Round off error of f(u) according to eq4.5: {roundOffError}")

plt.show()


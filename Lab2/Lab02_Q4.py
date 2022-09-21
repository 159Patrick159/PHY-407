#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def p(u):
    return (1-u)**8

def q(u):
    return 1 - 8*u + 28*u**2 - 56*u**3 + 70*u**4 - 56*u**5 + 28*u**6 - 8*u**7 + u**8

#%%
################################# Q4a #######################################
N = 500
u_range = np.linspace(0.98, 1.02, N)
p_vals = p(u_range)
q_vals = q(u_range)

plt.figure(figsize=(8,8))
plt.plot(u_range, p_vals, 'r.', label='p(u)')
plt.plot(u_range, q_vals, 'k.', label='q(u)')
plt.xlabel("u value", fontsize=16)
plt.ylabel("Polynomial numerical value", fontsize=16)
plt.title("Comparison of numerical value between expanded and simplified \nform of polynomial", fontsize=14)
plt.grid('on')
plt.legend()

#%%
################################# Q4b #######################################
diff_vals = p_vals - q_vals

plt.figure(figsize=(8,8))
plt.plot(u_range, diff_vals, ".", label='p(u)-q(u)')
plt.xlabel("u value", fontsize=16)
plt.ylabel("p(u)-q(u) numerical value", fontsize=16)
plt.title("Difference in numerical value between expanded and simplified \nform of polynomial", fontsize=14)
plt.grid('on')
plt.legend()

bin_array = np.linspace(min(diff_vals), max(diff_vals), 20)

plt.figure(figsize=(8,8), tight_layout=True)
plt.hist(diff_vals, bin_array, align='left', label='Data')
plt.xlabel("p(u)-q(u)",fontsize=16)
plt.ylabel("Number of occurance",fontsize=16)
plt.title("Distribution of numerical value of p(u)-q(u)",fontsize=16)

std_dev = np.std(diff_vals, ddof=1)

C = 1

term_p = (1-u_range)**8


mean_squared_diff = (sum(diff_vals**2))/N**2
error = C * np.sqrt(N) * np.sqrt(mean_squared_diff)

print(f"Standard deviation of distribution: {std_dev: .3e}")
print(f"Roundoff error estimate: {error: .3e}")

#%%
############################### Q4c #########################################
u_range = np.arange(0.9800, 0.9990, 0.0001)

p_vals = p(u_range)
q_vals = q(u_range)
err = abs(p_vals-q_vals)/abs(p_vals)

plt.figure(figsize=(8,8))
plt.plot(u_range, err, ".")
plt.xlabel("u value",fontsize=16)
plt.ylabel("Fractional Error?",fontsize=16)
plt.title("Fractional error", fontsize=16)
plt.yscale('log')
plt.grid('on')
#diff_vals = p_vals - q_vals
#std_dev = np.std(diff_vals, ddof=1)
#mean_squared_diff = (sum(u_range**2))/N
#error = C * np.sqrt(N) * np.sqrt(mean_squared_diff)
#print(f"Standard deviation of distribution: {std_dev: .3e}")
#print(f"Roundoff error estimate: {error: .3e}")



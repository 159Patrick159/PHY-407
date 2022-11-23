########################## 
# Contributors           # 
# Patrick Sandoval       #
# Kelvin Leong           #
##########################

################################# HEADER ##########################################
# This python script holds the code and plots for Q1 of Lab09 which solves the    #
# the time-dependent Schrodinger equation using Crank-Nicolson method             #
###################################################################################

# Import the necessary libary 
import numpy as np
import matplotlib.pyplot as plt
from MyFunctions import calc_normalization_constant, V, psi_init, calc_exp_X, calc_exp_E, calc_normalization_factor

# Constants
h_bar = 1.055e-34   #[m^2 kg/s]
L = 1e-8            #[m]
m = 9.109e-31       #[kg]
sigma = L/25.0
kappa = 500/L
x0 = L/5

# Spatial quantities
P = 1024
a = L/P
p_arr = np.arange(1, P, 1)
x_p_arr = p_arr*a - L/2

# Time quantities
N = 3000
tau = 1e-18
T = N*tau

# Define potential V(x)
potential_arr = V(x_p_arr, -L/2, L/2)

# Define the discretized Hamiltonian matrix
A = -h_bar**2/(2*m*a**2)
B = potential_arr - 2*A
Hamiltonian = np.zeros([len(p_arr), len(p_arr)])
for i in range(len(p_arr)):
    if i == 0:
        Hamiltonian[i][i] = B[i]
        Hamiltonian[i][i+1] = A
    elif i == len(p_arr)-1:
        Hamiltonian[i][i-1] = A
        Hamiltonian[i][i] = B[i]
    else:
        Hamiltonian[i][i-1] = A
        Hamiltonian[i][i] = B[i]
        Hamiltonian[i][i+1] = A

# Define the initial (t=0) wave fucntion from psi_1 to psi_{P-1}
# and normalize it (using relation in eq3)
psi_0 = psi_init(x_p_arr, x0, sigma, kappa)
nor_const = calc_normalization_constant(psi_init, P, -L/2, L/2, x0, sigma, kappa)
psi_0_nor = psi_0 / np.sqrt(nor_const)

# Construct L and R matrices used for Crank-Nicolson time-stepper method
L_matrix = np.eye(P-1) + 1j * tau/(2*h_bar) * Hamiltonian.copy()
L_matrix_inv = np.linalg.inv(L_matrix)
R_matrix = np.eye(P-1) - 1j * tau/(2*h_bar) * Hamiltonian.copy()

# Define the receiver of solution, and give initial condition
psi_all = np.zeros((N+1, P-1), dtype=complex)      #[Time][Spatial]
psi_all[0] = psi_0_nor

# Loop through time from n=0 to n=N-1 to solve for n+1 at each step
for n in range(N):
    psi_all[n+1, : ] = np.matmul(L_matrix_inv, np.matmul(R_matrix, psi_all[n, :]))

# End of solving the time-dependent Schrodinger equation
print("Finished solving time-dependent Schrodinger equation")
#%%
################################### Q1b ######################################
# Define position and time domain for plotting and calculating other quantities
x_domain = x_p_arr
t_domain = np.linspace(0, T, N+1)

# Define indices of times for plotting
n0, n1, n2, n3, n4 = 0, int(N/4), int(N/2), int(3*N/4), N

plt.figure(figsize=(10,10))
plt.plot(x_domain, np.real(psi_all[n0]), label=r'$\psi(t=0)$')
plt.plot(x_domain, np.real(psi_all[n1]), label=r'$\psi(t=T/4)$')
plt.plot(x_domain, np.real(psi_all[n2]), label=r'$\psi(t=T/2)$')
plt.plot(x_domain, np.real(psi_all[n3]), label=r'$\psi(t=3T/4)$')
plt.plot(x_domain, np.real(psi_all[n4]), label=r'$\psi(t=T)$')
plt.xlabel("Spatial Domain", fontsize=16)
plt.ylabel("Wave packet (real) amplitude", fontsize=16)
plt.title("The (real part) wavefunction at t=0, T/4, T/2, 3T/4, T", fontsize=16)
plt.grid('on')
plt.legend()
plt.savefig("Q1a.png")

# Calculate position expectation value
exp_X_t = calc_exp_X(psi_all, x_domain, a)
print("Finished calculating position expectation value")

plt.figure(figsize=(10, 10))
plt.plot(t_domain, exp_X_t)
plt.xlabel("Time domain", fontsize=16)
plt.ylabel("Expectation value <X>", fontsize=16)
plt.title("Expectation value of <X>(t)", fontsize=16)
plt.grid('on')
plt.savefig("Q1b.png")
#%%
################################### Q1c ######################################
# Calculate energy expectation value
exp_E_t = calc_exp_E(psi_all, Hamiltonian, a)
print("Finished calculating energy")

# Calculate normalization value over time
norm_t = calc_normalization_factor(psi_all, a)
print("Finished calculating normalization factor")

plt.figure(figsize=(10, 10))
plt.plot(t_domain, exp_E_t)
plt.xlabel("Time domain", fontsize=16)
plt.ylabel("Energy", fontsize=16)
plt.title("Energy evolution in the system", fontsize=16)
plt.ylim([0, 2e-17])
plt.grid('on')
plt.savefig("Q1c_1.png")

plt.figure(figsize=(10, 10))
plt.plot(t_domain, norm_t)
plt.xlabel("Time domain", fontsize=16)
plt.ylabel("Normalization factor", fontsize=16)
plt.title("Normalization of the wavefunction", fontsize=16)
plt.ylim([0.95, 1.05])
plt.grid('on')
plt.savefig("Q1c_2.png")
 
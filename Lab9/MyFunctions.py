#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 16:19:28 2022

@author: kelvinleong
"""
import numpy as np

def V(x_arr, a, b):
    mask = np.logical_or(x_arr < a, x_arr > b)
    mask = mask.astype(int) * 1e308     # Using largest floating point to represent infinity
    return mask
            

# Define Gaussian quad functions from textbook 
def gaussxw(N):
    # Initial approximation to roots of the Legendre polynomial
    a = np.linspace(3,4*N-1,N)/(4*N+2)
    x = np.cos(np.pi*a+1/(8*N*N*np.tan(a)))

    # Find roots using Newton's method
    epsilon = 1e-15
    delta = 1.0
    while delta>epsilon:
        p0 = np.ones(N,float)
        p1 = np.copy(x)
        for k in range(1,N):
            p0,p1 = p1,((2*k+1)*x*p1-k*p0)/(k+1)
        dp = (N+1)*(p0-x*p1)/(1-x*x)
        dx = p1/dp
        x -= dx
        delta = max(abs(dx))

    # Calculate the weights
    w = 2*(N+1)*(N+1)/(N*N*(1-x*x)*dp*dp)

    return x,w

def gaussxwab(N,a,b):
    x,w = gaussxw(N)
    return 0.5*(b-a)*x+0.5*(b+a),0.5*(b-a)*w

def calc_normalization_constant(f, N, a, b, x0, sigma, kappa):
    # Use gaussian quadrature to calculate the normalization integral
    # eq.3
    x, w = gaussxwab(N, a, b)
    total = 0
    for i in range(len(x)):
        total += w[i] * np.dot(np.conjugate(f(x[i],x0,sigma,kappa)), f(x[i],x0,sigma,kappa))
    return total

def psi(x, x0, sigma, kappa):
    return np.exp(-(x - x0)**2/(4*sigma**2) + kappa*x* 1j)


def calc_exp_X(psi, x, h):
    exp_X = np.zeros(len(psi))
    for n in range(len(psi)):
        exp_X[n] = simpsons(np.real(np.conj(psi[n]) * x * psi[n]), h)
        
    return exp_X

def calc_exp_E(psi, Hamiltonian, h):
    exp_E = np.zeros(len(psi))
    for n in range(len(psi)):
        integrand = np.matmul(np.conj(psi[n]), np.matmul(Hamiltonian, psi[n]))
        exp_E[n] = h * np.real(integrand) # Constant integration
        
    return exp_E

def calc_normalization_factor(psi, h):
    norm_t = np.zeros(len(psi))
    for n in range(len(psi)):
        norm_t[n] = simpsons(np.real(np.conj(psi[n]) * psi[n]), h)
        
    return norm_t

def simpsons(f, h):
    """
    Simpson integration method.
    """
    odd = 0
    even = 0
    for k in range(1, len(f) - 1):
        if k % 2 == 0:
            even += 4 * f[k]
        else:
            odd += 2 * f[k]

    return h / 3 * (f[0] + f[-1] + even + odd)


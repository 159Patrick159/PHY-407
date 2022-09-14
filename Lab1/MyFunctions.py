# -*- coding: utf-8 -*-
"""
Created on Wed Sep 14 10:15:32 2022

@Collaborator:  Kelvin Leong
                Patrick Sandoval
"""
# The file stores our own defined functions
import numpy as np

# Function to calculate the distance between Earth and Jupiter
def earth_jupiter_distance(x_E, y_E, x_J, y_J):
    """
    Parameters:
    x_E : [float] x-position of first object (Earth)
    y_E : [float] y-position of first object (Earth)
    x_J : [float] x-position of second object (Jupiter)
    y_J : [float] y-position of second object (Jupiter)
    
    Returns: [float] Euclidean distance between the two objects
    """
    dx = np.abs(x_E - x_J)
    dy = np.abs(y_E - y_J)
    return np.sqrt(dx**2 + dy**2)

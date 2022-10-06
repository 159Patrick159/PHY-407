import numpy as np
def make_wf(EigV,x):
    '''Given the eigenvector it calculates
    the corresponding wavefunction by the superposition
    principle'''
    LA = 5 # AA
    tmp = []
    for i in range(1,len(EigV)-1):
        tmp.append(EigV[i-1]*np.sin(i*np.pi*x/LA))
    return(sum(tmp))
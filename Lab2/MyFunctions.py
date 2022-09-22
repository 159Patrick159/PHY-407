def find_mean(data):
    '''Returns the mean of a data set'''
    import numpy as np
    return(np.sum(data)/len(data))
    
# Calculate area under the curve using trapezoidal rules (eq. 5.3)
def trapezoidal(f,a,b,N):
    # Width of slice
    h = (b - a) / N
    
    result = 0.5 * (f(a) + f(b))
    for k in range(1, N):
        result += f(a + k*h)
    result = h * result
    return result

# Calculate area under the curve using simpson's rules (eq. 5.9)
def simpsons(f,a,b,N):
    # Width of slice
    h = (b - a) / N
    
    odd_sum = 0
    even_sum = 0
    
    for k in range(1, N):
        if (k % 2)== 1:     # k is odd
            odd_sum += f(a + k*h)
        else:
            even_sum += f(a + k*h)
    result = h/3 * (f(a) + f(b) + 4 * odd_sum + 2 * even_sum)
    return result
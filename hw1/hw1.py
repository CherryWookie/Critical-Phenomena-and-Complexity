# Homework assignment for #4,5. This is the code for determining if a fixed point at given r,a,b is Unstable or Stable
# Michael Sell 2024

import numpy as np
import math

# Choose initial values r,a,b
r,a,b = 1,1,1

# Compute population values (N)
N1 = b + np.sqrt(r/a)
N2 = b - np.sqrt(r/a)

N_fp = [0,N2,N1]
fx_array = []
words = []

# Compute derivative of f(N') and use it to classify the fixed point
for N in N_fp:
    fx_prime = 3*(N**2) - 4*N*b + b**2 - (r/a)

    print(f'r,a,b: ', r,a,b)
    print(f'N:',N)
    print(f'd/dN: ', fx_prime)

    fx_array.append(fx_prime)
    if fx_prime < 0:
        print("UNSTABLE")
    elif fx_prime > 0:
        print("STABLE")
    else: 
        print("ZERO BOI")
    



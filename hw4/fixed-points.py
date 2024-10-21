import numpy as np
# Simple code to compute x,y pairs based on two equations for dx/dt, and dy/dt
# For Complexity homework #4
# Michael Sell September 2024

# Points vector
xyvals = [(1,-1),(1,1),(1,.02),(1,2),(0,0),(0,.5),(2,0),(.5,0),(0,2),(-1,0),(0,-1),(1,-3),(1,-.2)]
a = 2

# Compute x,y for given dx/dt, dy/dt values from vector
for x,y in xyvals:

    dxdt = (x*(1-x))
    dydt = y*(1-y) - a*x*y

    print(f"(x,y):({x},{y}) -> ({dxdt},{dydt})")
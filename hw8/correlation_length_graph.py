import numpy as np
import matplotlib.pyplot as plt

# Define parameters
p_c = 0.5927  # Critical percolation threshold
nu = 1        # Critical exponent for correlation length
A = 1         # Scaling constant for correlation length

# Define p values
p_left = np.linspace(0, p_c, 500, endpoint=False)  # Subcritical region
p_right = np.linspace(p_c, 1, 500, endpoint=False)  # Supercritical region

# Compute correlation length
xi_left = A * (p_c - p_left)**(-nu)
xi_right = A * (p_right - p_c)**(-nu)

print(f"Critical Value P_c: {p_c}")

# Plot the graph
plt.figure(figsize=(8, 6))
plt.plot(p_left, xi_left, label=r"$p < p_c$", color="blue")
plt.plot(p_right, xi_right, label=r"$p > p_c$", color="red")

# Add vertical line at p_c
plt.axvline(p_c, color="black", linestyle="--", label=r"$p = p_c$")

# Add labels and legend
plt.yscale("log")  # Correlation length diverges, use log scale
plt.title(r"Correlation Length $\xi$ as a Function of $p$", fontsize=14)
plt.xlabel(r"$p$ (Occupation Probability)", fontsize=12)
plt.ylabel(r"$\xi$ (Correlation Length)", fontsize=12)
plt.legend(fontsize=10)
plt.grid(alpha=0.3)

plt.show()

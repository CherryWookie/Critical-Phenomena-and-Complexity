import numpy as np
import matplotlib.pyplot as plt

# Define the renormalization function R_b(p)
def R_b(p):
    return p**4 - 4*p**3 + 4*p**2

# Define the diagonal line for fixed points
p_vals = np.linspace(0, 1, 500)
diagonal = p_vals

# Compute R_b(p)
R_vals = R_b(p_vals)

# Plot the graph
plt.figure(figsize=(8, 6))
plt.plot(p_vals, R_vals, label=r"$R_b(p) = p^4 - 4p^3 + 4p^2$", color="blue")
plt.plot(p_vals, diagonal, label=r"$R_b(p) = p$", color="red", linestyle="--")

# Highlight fixed points
fixed_points = [0, 1]  # trivial fixed points
plt.scatter(fixed_points, fixed_points, color="black", zorder=5, label="Trivial Fixed Points")

# Add labels and arrows to indicate flow direction
plt.arrow(0.2, 0.1, -0.1, -0.05, head_width=0.02, head_length=0.02, fc='green', ec='green')
plt.arrow(0.8, 0.9, 0.1, 0.05, head_width=0.02, head_length=0.02, fc='green', ec='green')
plt.text(0.15, 0.05, "Flow to $p=0$", fontsize=10, color="green")
plt.text(0.85, 0.95, "Flow to $p=1$", fontsize=10, color="green")

# Axes and labels
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.title("Renormalization Group Flow in $p$-Space", fontsize=14)
plt.xlabel(r"$p$ (Occupation Probability)", fontsize=12)
plt.ylabel(r"$R_b(p)$ (Renormalized Probability)", fontsize=12)
plt.legend(fontsize=10)
plt.grid(alpha=0.3)
plt.ylim(0, 1)
plt.xlim(0, 1)

plt.show()

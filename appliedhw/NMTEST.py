import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
from scipy.spatial import Delaunay

def generate_mesh(L, H, refinement=40):
    """Generate a triangular mesh for the rectangular domain."""
    x = np.linspace(0, L, refinement)
    y = np.linspace(0, H, refinement)
    xv, yv = np.meshgrid(x, y)
    points = np.c_[xv.ravel(), yv.ravel()]
    tri = Delaunay(points)
    return points, tri

def compute_element_matrix(vertices, k_over_mu):
    """Compute the local stiffness matrix for a triangular element."""
    x = vertices[:, 0]
    y = vertices[:, 1]
    A = 0.5 * np.abs(np.linalg.det(np.array([[1, x[0], y[0]], [1, x[1], y[1]], [1, x[2], y[2]]])))
    B = np.array([
        [x[1] - x[2], x[2] - x[0], x[0] - x[1]],
        [y[1] - y[2], y[2] - y[0], y[0] - y[1]]
    ]) / (2 * A)
    K = k_over_mu * A * (B.T @ B)
    return K

def assemble_global_matrix(points, tri, k_over_mu_function):
    """Assemble the global stiffness matrix."""
    n_points = len(points)
    K_global = lil_matrix((n_points, n_points))
    
    for simplex in tri.simplices:
        vertices = points[simplex]
        k_over_mu = k_over_mu_function(vertices.mean(axis=0))
        K_local = compute_element_matrix(vertices, k_over_mu)
        for i, ni in enumerate(simplex):
            for j, nj in enumerate(simplex):
                K_global[ni, nj] += K_local[i, j]
    
    return K_global

def apply_boundary_conditions(K, F, points, tri, vin, L, H, lin, lout):
    """Apply boundary conditions to the system."""
    n_points = len(points)
    for i, (x, y) in enumerate(points):
        if np.isclose(x, 0) and H-lin <= y <= H:  # Inflow boundary
            F[i] = -vin[0] * (1 / 1)  # Use vin and k/mu
        elif np.isclose(x, L) and 0 <= y <= lout:  # Outflow boundary
            K[i, :] = 0
            K[i, i] = 1
            F[i] = 0
        elif np.isclose(x, 0) or np.isclose(x, L) or np.isclose(y, 0) or np.isclose(y, H):  # Neumann (no flow)
            continue  # Natural Neumann condition

    return K, F

def solve_pressure(L, H, lin, lout, vin, k_over_mu_function, refinement=20):
    """Solve for the pressure field."""
    points, tri = generate_mesh(L, H, refinement)
    K_global = assemble_global_matrix(points, tri, k_over_mu_function)
    F_global = np.zeros(len(points))
    
    K_global, F_global = apply_boundary_conditions(K_global, F_global, points, tri, vin, L, H, lin, lout)
    
    # Solve the linear system
    P = spsolve(K_global.tocsr(), F_global)
    return points, tri, P

def compute_velocity(points, tri, P, k_over_mu_function):
    """Compute the velocity field v = -k/(mu*epsilon) * grad(p)."""
    U = np.zeros(len(tri.simplices))
    V = np.zeros(len(tri.simplices))
    point_tri = np.zeros((len(tri.simplices), 2))

    for l, simplex in enumerate(tri.simplices):
        j1, j2, j3 = simplex
        x1, y1, f1 = points[j1, 0], points[j1, 1], P[j1]
        x2, y2, f2 = points[j2, 0], points[j2, 1], P[j2]
        x3, y3, f3 = points[j3, 0], points[j3, 1], P[j3]

        A = ((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)) / 2
        dphi1 = 1 / (2 * A) * np.array([y2 - y3, x3 - x2])
        dphi2 = 1 / (2 * A) * np.array([y3 - y1, x1 - x3])
        dphi3 = 1 / (2 * A) * np.array([y1 - y2, x2 - x1])

        el_vel = f1 * dphi1 + f2 * dphi2 + f3 * dphi3
        k_over_mu = k_over_mu_function(np.mean(points[simplex], axis=0))

        U[l] = k_over_mu * el_vel[0]
        V[l] = k_over_mu * el_vel[1]

        point_tri[l, 0] = (x1 + x2 + x3) / 3
        point_tri[l, 1] = (y1 + y2 + y3) / 3

    return point_tri, U, V

def k_over_mu_function(point):
    x, y = point
    if (x - circle_center[0])**2 + (y - circle_center[1])**2 < circle_radius**2:  # Inside the circular region
        return value_inside_circle
    else:
        return 1
    
# Define the parameters
L, H = 2, 1
lin, lout = 0.1, 0.2
vin = (1, 0)
circle_center=(L/2, H/2)
circle_radius=0.2
value_inside_circle = 0.1

# Solve the problem
points, tri, P = solve_pressure(L, H, lin, lout, vin, k_over_mu_function)

# Plot the mesh
#plot_mesh(points, tri)

plt.figure(1,figsize=(8,5))
plt.tricontourf(points[:, 0], points[:, 1], tri.simplices, P, levels=50)
plt.colorbar(label="Pressure p")

# Plot the domain boundaries
plt.plot([0, 0, L, L, 0], [0, H, H, 0, 0], color='black', linestyle='--', label='Domain Boundary')
plt.plot([0, 0], [H - lin,H], color='green', linewidth=4, label='Inflow Boundary')
plt.plot([L, L], [0, lout], color='red', linewidth=4, label='Outflow Boundary')

# Plot the circular region
circle = plt.Circle(circle_center, circle_radius, color='blue', linestyle='--', fill=False, linewidth=1, label='Low k/μ Region')
plt.gca().add_artist(circle)

plt.title("Pressure Field")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
#plt.legend()
plt.gca().set_aspect('equal')
plt.show(block=False)

point_tri, U,V = compute_velocity(points, tri, P, k_over_mu_function)

# Plot velocity field
plt.figure(2,figsize=(8,5))
plt.quiver(point_tri[:, 0], point_tri[:, 1], U, V, color='blue', scale=40, alpha=0.6)
plt.title("Velocity Field")
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.gca().set_aspect('equal')
plt.show()

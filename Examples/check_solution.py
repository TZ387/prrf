import numpy as np

"""
This Python code solves the problem of determining the electric potential distribution 
between two infinite plates that are separated by a distance L = 10.5 mm, with four 
semi-infinite regions of different conductivities between them. The setup is as follows:

- One plate is at potential V = 0 V (z = 0) and the other plate is at potential V = 10 V (z = L).
- Between the plates, there are four semi-infinite regions with different conductivities (σ) and lengths (L):
    - Region 1: σ1 = 0.3 S/m, length L1 = 0.5 mm
    - Region 2: σ2 = 0.22 S/m, length L2 = 1.5 mm
    - Region 3: σ3 = 0.025 S/m, length L3 = 3.5 mm
    - Region 4: σ4 = 0.5 S/m, length L4 = 5.0 mm

The goal is to solve for the potential in each of these regions by assuming the potential in each region is a linear function 
of the form V_i(z) = k_i * z + C_i, where k_i and C_i are the slope and intercept specific to region i.

To solve for the constants (k_i, C_i):
- We enforce continuity of potential and current at the interfaces between the regions (z = L1, z = L1 + L2, etc.).
- The current continuity at the interfaces means that the current flowing from one region to the next must be the same, 
which translates into the relation σ_i * k_i = σ_(i+1) * k_(i+1) for the slopes of the potential in adjacent regions.
- Additionally, the boundary conditions are that the potential is 0 V at z = 0 and 10 V at z = L.

The system of equations is solved using NumPy's linear algebra solver. After solving, the code calculates the potential 
at the edges between the regions and prints the potentials at the interface points: z = L1, z = L1 + L2, and z = L1 + L2 + L3.
"""


# Given data
sigma1 = 0.3   # S/m
sigma2 = 0.22  # S/m
sigma3 = 0.025 # S/m
sigma4 = 0.5   # S/m

L1 = 0.5e-3    # m
L2 = 1.5e-3    # m
L3 = 3.5e-3    # m
L4 = 5.0e-3    # m
L = L1 + L2 + L3 + L4  # Total length

# Potentials at the boundaries
V0 = 0        # V at z = 0
V4 = 10       # V at z = L

# System of equations for the potentials
# We need to solve for k1, k2, k3, k4, and C2, C3, C4 (since C1 = 0, V0 is at z=0)

# Current continuity relations:
# sigma1 * k1 = sigma2 * k2
# sigma2 * k2 = sigma3 * k3
# sigma3 * k3 = sigma4 * k4

# Create a system of linear equations:
A = np.array([
    [L1, -L1, 0, 0, -1, 0, 0],                # Continuity of potential at z = L1
    [0, L1 + L2, -L1 - L2, 0, 1, -1, 0],           # Continuity of potential at z = L1 + L2
    [0, 0, L1 + L2 + L3, -L1 - L2 - L3, 0, 1, -1],      # Continuity of potential at z = L1 + L2 + L3
    [sigma1, -sigma2, 0, 0, 0, 0, 0],       # Continuity of current at z = L1
    [0, sigma2, -sigma3, 0, 0, 0, 0],       # Continuity of current at z = L1 + L2
    [0, 0, sigma3, -sigma4, 0, 0, 0],       # Continuity of current at z = L1 + L2 + L3
    [0, 0, 0, L, 0, 0, 1]                   # Boundary condition V(L) = 10 V
])

# The right-hand side (known potentials and boundary conditions)
b = np.array([0, 0, 0, 0, 0, 0, V4])

# Solve the system
solution = np.linalg.solve(A, b)

# Extract the potentials and slopes
k1, k2, k3, k4 = solution[0:4]
C2, C3, C4 = solution[4:7]

# Calculate the potentials at the edges between plates
V_L1 = k1 * L1           # Potential at z = L1
V_L1_L2 = k2 * (L1 + L2) + C2  # Potential at z = L1 + L2
V_L1_L2_L3 = k3 * (L1 + L2 + L3) + C3  # Potential at z = L1 + L2 + L3

# Print the results
print(f"Potential at z = L1: {V_L1:.4f} V")
print(f"Potential at z = L1 + L2: {V_L1_L2:.4f} V")
print(f"Potential at z = L1 + L2 + L3: {V_L1_L2_L3:.4f} V")
print(f"Potential at z = L (boundary): {V4} V")

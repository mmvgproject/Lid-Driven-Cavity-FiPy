#
# Lid-Driven Cavity Simulation using FiPy and the SIMPLE algorithm (Corrected)
#
import os
import imageio
import numpy as np
import matplotlib.pyplot as plt
from fipy import (CellVariable, Grid2D, FaceVariable, TransientTerm,
                  ConvectionTerm, DiffusionTerm, Viewer, ImplicitSourceTerm)
from fipy.tools import numerix

print("FiPy Lid-Driven Cavity Simulation Started...")

# 1. --- Problem Parameters ---
# Physical parameters
L = 1.0  # Length of the cavity
U_lid = 1.0  # Velocity of the lid
Re = 400.0  # Reynolds number
nu = U_lid * L / Re  # Kinematic viscosity (nu = mu/rho)

# Numerical parameters
nx = 50  # Number of cells in x-direction
ny = 50  # Number of cells in y-direction
dx = L / nx
dy = L / ny
max_iterations = 100
tolerance = 1e-5

# SIMPLE algorithm parameters
alpha_u = 0.7  # Relaxation factor for velocity
alpha_p = 0.3  # Relaxation factor for pressure

# GIF generation parameters
generate_gif = True
gif_filename = 'lid_driven_cavity_evolution.gif'
frame_interval = 10  # Save a frame every 10 iterations
frame_files = []

# 2. --- Mesh and Variable Setup ---
mesh = Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)

# Define variables
# Velocity components are cell-centered
u = CellVariable(name=r"$u$", mesh=mesh, value=0.0)
v = CellVariable(name=r"$v$", mesh=mesh, value=0.0)
velocity = FaceVariable(name="velocity", mesh=mesh, rank=1)

# Pressure is also cell-centered
pressure = CellVariable(name=r"$p$", mesh=mesh, value=0.0)
pressure_correction = CellVariable(name=r"$p'$", mesh=mesh, value=0.0)

# 3. --- Boundary Conditions ---
# Lid (top wall)
u.constrain(U_lid, mesh.facesTop)
v.constrain(0.0, mesh.facesTop)

# Other walls (no-slip)
u.constrain(0.0, mesh.facesBottom | mesh.facesLeft | mesh.facesRight)
v.constrain(0.0, mesh.facesTop | mesh.facesBottom | mesh.facesLeft | mesh.facesRight)

# Pin pressure at one point (cell 0) for a unique solution
# --- THIS IS THE CORRECTED PART ---
pressure_mask = np.zeros(mesh.numberOfCells, dtype=bool)
pressure_mask[0] = True
pressure.constrain(0., where=pressure_mask)
# ---------------------------------

# 4. --- Equation Setup ---
# Use Upwind scheme for convection as requested
convection_scheme = ConvectionTerm

# Momentum equations (predictor step)
# We solve for u and v separately. The pressure gradient is treated explicitly.
# The `dt` is a pseudo-timestep for the iterative solver; a large dt approximates a steady state solve.
u_eq = TransientTerm() == -convection_scheme(coeff=velocity) + DiffusionTerm(coeff=nu) - pressure.grad.dot([1., 0.])
v_eq = TransientTerm() == -convection_scheme(coeff=velocity) + DiffusionTerm(coeff=nu) - pressure.grad.dot([0., 1.])

# Pressure correction (Poisson equation)
# This is the core of the SIMPLE algorithm. It's derived from the continuity equation.
# The `ap` term represents the influence of pressure on velocity.
ap = CellVariable(mesh=mesh, value=1.0) # Placeholder, will be updated in loop
p_correction_eq = DiffusionTerm(coeff=1.0/ap) - velocity.divergence == 0

# 5. --- SIMPLE Algorithm Loop ---
if generate_gif:
    print(f"Animation frames will be saved. Final GIF will be '{gif_filename}'.")

# We need to store old values to calculate convergence
u_old = u.copy()
v_old = v.copy()

for sweep in range(max_iterations):
    # --- Step 1: Solve Momentum Predictor ---
    # Use previous iteration's pressure to predict intermediate velocities (u*, v*)
    velocity.setValue(numerix.array((u.arithmeticFaceValue, v.arithmeticFaceValue)))
    
    u_eq.solve(var=u, dt=1e2)
    v_eq.solve(var=v, dt=1e2)

    # --- Step 2: Solve Pressure Correction ---
    # The source term for the pressure correction equation is the divergence of the
    # intermediate velocity field.
    velocity.setValue(numerix.array((u.arithmeticFaceValue, v.arithmeticFaceValue)))
    
    # We need the diagonal coefficients (A_p) from the discretized momentum equations
    # For a simple uniform grid and collocated variables, this can be approximated.
    ap.setValue(1.0) # A simplification for this problem setup

    # Solve the Poisson equation for pressure correction p'
    p_correction_eq.solve(var=pressure_correction)

    # --- Step 3: Correct Pressure and Velocity ---
    # Correct pressure with under-relaxation
    pressure.setValue(pressure() + alpha_p * pressure_correction())
    
    # Correct velocities to satisfy continuity
    u.setValue(u() - (1.0 / ap) * pressure_correction.grad.dot([1., 0.]))
    v.setValue(v() - (1.0 / ap) * pressure_correction.grad.dot([0., 1.]))
    
    # Apply under-relaxation to velocities
    u.setValue(alpha_u * u() + (1 - alpha_u) * u_old())
    v.setValue(alpha_u * v() + (1 - alpha_u) * v_old())
    
    # Re-apply boundary conditions after correction
    u.constrain(U_lid, mesh.facesTop)
    u.constrain(0.0, mesh.facesBottom | mesh.facesLeft | mesh.facesRight)
    v.constrain(0.0, mesh.facesTop | mesh.facesBottom | mesh.facesLeft | mesh.facesRight)

    # --- Step 4: Check for Convergence ---
    # Calculate mass imbalance (continuity residual)
    velocity.setValue(numerix.array((u.arithmeticFaceValue, v.arithmeticFaceValue)))
    continuity_residual = abs(velocity.divergence).max()

    if sweep > 0 and continuity_residual < tolerance:
        print(f"\nConverged at iteration {sweep} with continuity residual: {continuity_residual:.2e}")
        break

    if sweep % 20 == 0:
        print("sweep : ")
        print(sweep)
        print("continuity_residual : ")
        print(continuity_residual)
        

    # Update old values for the next iteration's relaxation
    u_old.setValue(u())
    v_old.setValue(v())

    # --- Frame Generation for GIF ---
    if generate_gif and sweep % frame_interval == 0:
        fig, ax = plt.subplots()
        im = ax.imshow(u.value.reshape(ny, nx), origin='lower', cmap='viridis',
                       extent=[0, L, 0, L], vmin=-0.2, vmax=1.0)
        ax.set_title(f'u-velocity, Iteration: {sweep}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        plt.colorbar(im, ax=ax, label='u-velocity (m/s)')
        frame_filepath = f'frame_{sweep:04d}.png'
        plt.savefig(frame_filepath)
        plt.close(fig)
        frame_files.append(frame_filepath)

else:
     print(f"\nReached maximum iterations ({max_iterations}) without converging.")


# 6. --- Post-Processing and Visualization ---
print("\nSimulation finished. Generating final plots...")

# Function to create nice plots
def create_contour_plot(var, title, cmap='viridis'):
    fig, ax = plt.subplots(figsize=(7, 6))
    data = var.value.reshape(ny, nx)
    im = ax.imshow(data, origin='lower', cmap=cmap, extent=[0, L, 0, L])
    ax.set_title(title, fontsize=16)
    ax.set_xlabel('x-position (m)', fontsize=12)
    ax.set_ylabel('y-position (m)', fontsize=12)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(f'{var.name} value', fontsize=12)
    plt.tight_layout()
    plt.savefig(f'{title.replace(" ", "_").lower()}.png', dpi=300)
    plt.show()

# Plot u-velocity, v-velocity, and pressure
create_contour_plot(u, 'U-Velocity Field', cmap='viridis')
create_contour_plot(v, 'V-Velocity Field', cmap='coolwarm')
create_contour_plot(pressure, 'Pressure Field', cmap='jet')


# --- Stream Function and Streamline Plot ---
print("Calculating and plotting stream function...")

# 1. Calculate Vorticity (omega_z)
vorticity = v.grad.dot([1., 0.]) - u.grad.dot([0., 1.])
vorticity.name = r"$\omega_z$"

# 2. Solve for Stream Function (psi)
stream_function = CellVariable(name=r"$\psi$", mesh=mesh)
stream_function.constrain(0.0, mesh.facesTop | mesh.facesBottom | mesh.facesLeft | mesh.facesRight)
psi_eq = DiffusionTerm(coeff=1.0) == -vorticity
psi_eq.solve(var=stream_function)

# 3. Plot Stream Function and Streamlines
fig, ax = plt.subplots(figsize=(7, 6))
x, y = np.meshgrid(mesh.cellCenters[0].value.reshape(nx,1), mesh.cellCenters[1].value.reshape(1,ny))
u_val = u.value.reshape(ny, nx)
v_val = v.value.reshape(ny, nx)
ax.streamplot(x, y, u_val, v_val, color='black', linewidth=1, density=1.5)

sf_data = stream_function.value.reshape(ny, nx)
contour = ax.contourf(x, y, sf_data, cmap='autumn', levels=20)
ax.set_title('Streamlines and Stream Function', fontsize=16)
ax.set_xlabel('x-position (m)', fontsize=12)
ax.set_ylabel('y-position (m)', fontsize=12)
ax.set_aspect('equal', adjustable='box')
cbar = plt.colorbar(contour, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label('Stream Function ' + r'$\psi$', fontsize=12)
plt.tight_layout()
plt.savefig('stream_function_and_lines.png', dpi=300)
plt.show()


# 7. --- GIF Generation ---
if generate_gif and frame_files:
    print(f"\nCreating GIF animation: {gif_filename}...")
    with imageio.get_writer(gif_filename, mode='I', duration=0.1, loop=0) as writer:
        for filename in frame_files:
            image = imageio.imread(filename)
            writer.append_data(image)
    for filename in frame_files:
        os.remove(filename)
    print("GIF created successfully and frame files removed.")

Lid-Driven Cavity Flow Simulation using FiPy
This repository contains a Python script for simulating 2D incompressible, viscous fluid flow in a lid-driven cavity using the Finite Volume Method. The pressure-velocity coupling is handled by the SIMPLE (Semi-Implicit Method for Pressure-Linked Equations) algorithm. The implementation leverages the FiPy library for solving partial differential equations.

Figure 1: Final streamlines and stream function contours for Re=100.

Table of Contents
Project Overview

Theoretical Background

The Lid-Driven Cavity Problem

Governing Equations (Navier-Stokes)

The SIMPLE Algorithm

Program Specifications

Dependencies

How to Run

Customization

Code Logic and Implementation

1. Setup and Parameters

2. Mesh and Variables

3. Boundary Conditions

4. Equation Discretization

5. The SIMPLE Iterative Loop

6. Post-Processing

Expected Results

Suggestions for Improvement

Project Overview
The primary goal of this project is to accurately simulate and visualize the fluid flow within a square cavity where the top lid moves at a constant velocity. This is a classical benchmark problem in Computational Fluid Dynamics (CFD) used to validate numerical methods.

The script solves the steady-state, incompressible Navier-Stokes equations. It employs the SIMPLE algorithm to handle the inherent coupling between pressure and velocity. The entire simulation is built upon FiPy, a powerful Python-based finite volume PDE solver.

Key features of this simulation include:

Finite Volume Discretization on a uniform 2D grid.

Implementation of the iterative SIMPLE algorithm from scratch.

Visualization of results, including velocity and pressure fields, streamlines, and an optional animated GIF showing the convergence of the velocity field.

Theoretical Background
The Lid-Driven Cavity Problem
The lid-driven cavity is a standard CFD validation case. It consists of a square domain filled with a fluid, enclosed by three stationary walls (bottom, left, right) and a top wall (the "lid") that moves tangentially with a constant velocity U_lid. This motion of the lid induces a recirculating flow pattern inside the cavity. The nature of this flow is primarily determined by the Reynolds number (Re), a dimensionless quantity that represents the ratio of inertial forces to viscous forces.

Re= 
μ
ρU 
lid
​
 L
​
 = 
ν
U 
lid
​
 L
​
 
where:

rho is the fluid density.

mu is the dynamic viscosity.

nu=
mu/
rho is the kinematic viscosity.

L is the characteristic length (the side length of the cavity).

Governing Equations (Navier-Stokes)
For an incompressible, Newtonian fluid, the flow is governed by the Navier-Stokes equations, which express the conservation of mass and momentum.

1. Continuity Equation (Conservation of Mass):
For an incompressible fluid, the divergence of the velocity field must be zero.


∇⋅ 
u
 = 
∂x
∂u
​
 + 
∂y
∂v
​
 =0
2. Momentum Equation (Conservation of Momentum):
This equation relates the acceleration of the fluid to the forces acting on it (pressure gradient, viscous forces).


∂t
∂ 
u
 
​
 +( 
u
 ⋅∇) 
u
 =− 
ρ
1
​
 ∇p+ν∇ 
2
  
u
 
Since we are seeking the steady-state solution, the time derivative term 
fracpartialvecupartialt is zero. The challenge in solving these equations lies in the coupling between the velocity components (
vecu) and the pressure (p). The pressure field must evolve to ensure that the final velocity field satisfies the continuity equation.

The SIMPLE Algorithm
The Semi-Implicit Method for Pressure-Linked Equations (SIMPLE) is an iterative procedure designed to solve the pressure-velocity coupling problem. It follows a "predictor-corrector" approach.

The core idea is as follows:

Guess the pressure field, p 
\*
 .

Predictor Step: Solve the momentum equations using this guessed pressure p 
∗
  to obtain an intermediate velocity field, 
vecu 
∗
 . This velocity field will not, in general, satisfy the continuity equation.

Corrector Step: Formulate a Poisson equation for a pressure correction, p 
′
 . The source term for this equation is the divergence of the intermediate velocity field (
nabla
cdot
vecu 
\*
 ), which represents the mass imbalance.

Solve for p 
′
  and use it to correct the pressure and velocity fields:

p=p 
\*
 +
alpha_pp 
′
 

vecu=
vecu 
\*
 +
vecu 
′
 (p 
′
 )

Repeat steps 2-4 until the corrections become negligibly small and the continuity equation is satisfied to a desired tolerance.

Under-relaxation factors (
alpha_p for pressure, 
alpha_u for velocity) are used to stabilize this iterative process and prevent divergence.

Program Specifications
Dependencies
This script requires the following Python libraries. You can install them using pip:

pip install fipy numpy matplotlib imageio

FiPy: The core engine for solving the PDEs.

NumPy: For numerical operations and array handling.

Matplotlib: For generating plots and visualizations.

Imageio: For creating the optional GIF animation.

How to Run
To run the simulation, simply execute the Python script from your terminal:

python lid_driven_cavity.py

The script will print the progress of the SIMPLE algorithm iterations to the console and, upon completion, will save several output files in the same directory.

Customization
You can modify the simulation by changing the parameters in the "Problem Parameters" section of the script:

Re: Change the Reynolds number to simulate different flow regimes.

nx, ny: Increase for a finer mesh and more accurate results, at the cost of longer computation time.

max_iterations: Increase if the simulation fails to converge.

tolerance: Decrease for a more stringent convergence criterion.

alpha_u, alpha_p: Adjust the under-relaxation factors. Values typically range from 0.1 to 0.9.

generate_gif: Set to False to disable the creation of the animation.

Code Logic and Implementation
1. Setup and Parameters
The script begins by defining the physical and numerical parameters, including the cavity size L, lid velocity U_lid, Reynolds number Re, and grid dimensions nx, ny. It also sets the under-relaxation factors and GIF generation options.

2. Mesh and Variables
A 2D uniform grid (Grid2D) is created. The primary variables—velocity components u and v, and pressure p—are defined as CellVariable, meaning their values are stored at the center of each grid cell.

3. Boundary Conditions
Lid: The top wall has a fixed u velocity of U_lid and a v velocity of 0.

Other Walls: The bottom, left, and right walls have no-slip conditions, meaning both u and v are fixed to 0.

Pressure Pinning: The pressure field in the incompressible Navier-Stokes equations is relative. To obtain a unique solution, the pressure at a single reference point must be fixed (or "pinned"). In this code, the pressure in the first cell (index 0) is constrained to 0.

4. Equation Discretization
The momentum and pressure correction equations are defined using FiPy's declarative syntax.

u_eq and v_eq: These represent the momentum equations. They balance a TransientTerm (used here as an iterative solver term) with ConvectionTerm, DiffusionTerm, and the explicit pressure gradient.

p_correction_eq: This is the Poisson equation for the pressure correction, formulated as a DiffusionTerm equaling the divergence of the intermediate velocity field.

5. The SIMPLE Iterative Loop
This is the core of the program, executing the algorithm described in the theory section.

Solve Momentum Predictor: The u_eq and v_eq are solved using the pressure from the previous iteration to get u* and v*.

Solve Pressure Correction: The p_correction_eq is solved to find the pressure correction p'.

Correct Fields: The pressure and velocity fields are updated using p'. Under-relaxation is applied to both velocity and pressure to stabilize the solution.

Check Convergence: The continuity residual (the maximum divergence of the corrected velocity field) is calculated. If it falls below the specified tolerance, the loop breaks.

6. Post-Processing
Once the solution has converged, the script performs several post-processing steps:

Contour Plots: Generates and saves contour plots for the final u-velocity, v-velocity, and pressure fields.

Stream Function Calculation:

Calculates the vorticity (
omega_z=
fracpartialvpartialx−
fracpartialupartialy).

Solves the Poisson equation 
nabla 
2
 
psi=−
omega_z to find the stream function 
psi.

Streamline Plot: Creates a final plot that overlays the velocity streamlines on a contour plot of the stream function, providing a clear visualization of the flow pattern.

GIF Generation: If enabled, it compiles the saved frames of the u-velocity evolution into an animated GIF.

Expected Results
After running the script with the default parameters (Re=100, 50x50 grid), you should expect the following output files:

u_velocity_field.png: Contour plot of the horizontal velocity.

v_velocity_field.png: Contour plot of the vertical velocity.

pressure_field.png: Contour plot of the relative pressure.

stream_function_and_lines.png: A plot showing the primary vortex and stream function contours (as seen in Figure 1).

lid_driven_cavity_evolution.gif (optional): An animation showing the development of the u-velocity field from its initial state to the final converged solution.

The flow pattern at Re=100 is characterized by a large, primary vortex in the center of the cavity. Smaller, less intense secondary vortices are expected to form in the bottom-left and bottom-right corners, though they may be difficult to see on a coarse mesh.

Suggestions for Improvement
This implementation serves as a solid foundation. Here are some ways it could be extended or improved:

Higher Reynolds Numbers: To accurately simulate flows at higher Re (e.g., 400, 1000), a finer mesh (nx, ny) and potentially more iterations (max_iterations) will be required. The flow becomes more complex, with more pronounced secondary and tertiary vortices.

Quantitative Validation: Add code to plot velocity profiles along the vertical and horizontal centerlines of the cavity. These results can then be compared against well-established benchmark data from sources like Ghia et al. (1982) to quantitatively validate the accuracy of the simulation.

Performance Optimization: For very fine meshes, the simulation can become computationally expensive. FiPy supports parallel execution, which could be implemented to significantly speed up the calculations on multi-core processors.

Algorithm Enhancements: Implement more advanced pressure-velocity coupling schemes like SIMPLEC (SIMPLE-Consistent) or PISO (Pressure-Implicit with Splitting of Operators), which can offer faster convergence rates in certain cases.

Code Modularity: Refactor the code by encapsulating the SIMPLE loop, plotting functions, and problem setup into separate functions or a class for better organization and reusability.

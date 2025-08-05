# Lid-Driven Cavity Simulation using FiPy and SIMPLE Algorithm

## Project Title
**Lid-Driven Cavity Simulation with FiPy and SIMPLE Algorithm**

This project implements a numerical simulation of the lid-driven cavity flow, a classic benchmark problem in computational fluid dynamics (CFD). The simulation uses the **FiPy** library, a Python-based finite volume PDE solver, and the **SIMPLE (Semi-Implicit Method for Pressure-Linked Equations)** algorithm to solve the incompressible Navier-Stokes equations. The code models a 2D square cavity with a moving lid, generating detailed velocity and pressure fields, vorticity, streamlines, and an animated GIF of the solution evolution.

---

## Project and Program Specifications

### Overview
The lid-driven cavity problem simulates the flow of an incompressible fluid in a square domain where the top boundary (the "lid") moves at a constant velocity, while the other three walls are stationary. This creates a complex flow pattern driven by viscous effects and pressure gradients, making it an ideal test case for CFD algorithms.

### Program Specifications
- **Programming Language**: Python 3.x
- **Key Libraries**:
  - **FiPy**: For solving partial differential equations (PDEs) using the finite volume method.
  - **NumPy**: For numerical operations and array handling.
  - **Matplotlib**: For visualization of velocity, pressure, and streamline plots.
  - **ImageIO**: For generating animated GIFs of the solution evolution.
- **Algorithm**: SIMPLE (Semi-Implicit Method for Pressure-Linked Equations).
- **Numerical Method**: Finite volume method with an upwind convection scheme.
- **Boundary Conditions**:
  - Top wall (lid): Constant velocity \( u = U_{\text{lid}}, v = 0 \).
  - Other walls (left, right, bottom): No-slip condition (\( u = 0, v = 0 \)).
  - Pressure pinned at one point to ensure a unique solution.
- **Physical Parameters**:
  - Cavity length: \( L = 1.0 \, \text{m} \).
  - Lid velocity: \( U_{\text{lid}} = 1.0 \, \text{m/s} \).
  - Reynolds number: \( \text{Re} = 100 \).
  - Kinematic viscosity: \( \nu = \frac{U_{\text{lid}} \cdot L}{\text{Re}} \).
- **Numerical Parameters**:
  - Grid size: \( 50 \times 50 \) cells.
  - Relaxation factors: \( \alpha_u = 0.7 \) (velocity), \( \alpha_p = 0.3 \) (pressure).
  - Convergence tolerance: \( 10^{-5} \).
  - Maximum iterations: 100.
- **Output**:
  - Contour plots for \( u \)-velocity, \( v \)-velocity, and pressure fields.
  - Streamline plot with stream function contours.
  - Optional animated GIF showing the evolution of the \( u \)-velocity field.

### Dependencies
To run the code, install the required Python packages:
```bash
pip install fipy numpy matplotlib imageio
```

---

## Project Goals
The primary goals of this project are:
1. **Accurate Simulation**: Solve the 2D incompressible Navier-Stokes equations for the lid-driven cavity problem using the SIMPLE algorithm.
2. **Visualization**: Generate high-quality visualizations of the velocity, pressure, and stream function fields to analyze the flow behavior.
3. **Educational Tool**: Provide a well-documented codebase to demonstrate the implementation of the SIMPLE algorithm in FiPy.
4. **Animation**: Create an animated GIF to visualize the iterative convergence of the solution.
5. **Reproducibility**: Ensure the code is modular, well-commented, and suitable for use in academic and research settings.

---

## Logic Used in the Project

### Overview of the SIMPLE Algorithm
The SIMPLE algorithm is used to couple the velocity and pressure fields in the incompressible Navier-Stokes equations. It iteratively solves the momentum equations and a pressure correction equation to ensure continuity (zero divergence of the velocity field). The key steps are:

1. **Momentum Predictor**: Solve the momentum equations using the current pressure field to obtain intermediate velocities (\( u^*, v^* \)).
2. **Pressure Correction**: Solve a Poisson equation for the pressure correction (\( p' \)) to enforce continuity.
3. **Velocity and Pressure Update**: Correct the velocities and pressure using under-relaxation to stabilize convergence.
4. **Convergence Check**: Monitor the continuity residual (velocity divergence) and stop when it falls below a specified tolerance.

### Code Structure
The code is organized into seven main sections:
1. **Problem Parameters**: Defines physical and numerical parameters (e.g., Reynolds number, grid size, relaxation factors).
2. **Mesh and Variable Setup**: Creates a 2D grid using FiPy and initializes cell-centered variables for velocity (\( u, v \)), pressure (\( p \)), and pressure correction (\( p' \)).
3. **Boundary Conditions**: Applies the lid velocity and no-slip conditions, with pressure pinned at one point.
4. **Equation Setup**: Defines the momentum equations (with convection and diffusion terms) and the pressure correction equation.
5. **SIMPLE Algorithm Loop**: Iteratively solves the momentum and pressure correction equations, applies under-relaxation, and checks for convergence.
6. **Post-Processing and Visualization**: Generates contour plots for velocity, pressure, and stream function, along with streamlines.
7. **GIF Generation**: Creates an animated GIF of the \( u \)-velocity field evolution.

### Key Implementation Details
- **Finite Volume Method**: The FiPy library discretizes the PDEs using a finite volume approach, ensuring conservation of mass and momentum.
- **Upwind Convection Scheme**: Used for the convection term to improve numerical stability.
- **Under-Relaxation**: Applied to velocity and pressure updates to prevent oscillations and ensure convergence.
- **Pressure Pinning**: Fixes the pressure at one cell to eliminate the null space in the pressure Poisson equation.
- **Stream Function**: Calculated by solving a Poisson equation for the stream function (\( \psi \)) using vorticity as the source term.

---

## Theory Used in the Project

### Governing Equations
The lid-driven cavity flow is governed by the 2D incompressible Navier-Stokes equations:
1. **Continuity Equation** (mass conservation):
   \[
   \nabla \cdot \mathbf{u} = 0
   \]
   where \( \mathbf{u} = (u, v) \) is the velocity vector.
2. **Momentum Equations** (momentum conservation):
   \[
   \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla) \mathbf{u} = -\frac{1}{\rho} \nabla p + \nu \nabla^2 \mathbf{u}
   \]
   For steady-state, the time derivative is approximated using a pseudo-timestep in the SIMPLE algorithm.

### SIMPLE Algorithm Theory
The SIMPLE algorithm addresses the challenge of coupling velocity and pressure in incompressible flows. Since the continuity equation does not directly involve pressure, SIMPLE introduces a pressure correction to enforce \( \nabla \cdot \mathbf{u} = 0 \). The algorithm:
- Solves the momentum equations with an estimated pressure field.
- Computes the velocity divergence to form a Poisson equation for the pressure correction.
- Updates the velocity and pressure fields iteratively until convergence.

### Stream Function and Vorticity
The stream function (\( \psi \)) is introduced to visualize the flow field, defined such that:
\[
u = \frac{\partial \psi}{\partial y}, \quad v = -\frac{\partial \psi}{\partial x}
\]
The vorticity (\( \omega_z \)) is calculated as:
\[
\omega_z = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y}
\]
The stream function satisfies the Poisson equation:
\[
\nabla^2 \psi = -\omega_z
\]
This is solved to obtain \( \psi \), which is used to plot streamlines and contours.

### Reynolds Number
The Reynolds number (\( \text{Re} = \frac{U_{\text{lid}} L}{\nu} \)) governs the flow regime. At \( \text{Re} = 100 \), the flow is laminar with a single primary vortex and weak secondary vortices, making it suitable for this simulation.

---

## Expected Results
The simulation produces the following outputs:
1. **Velocity Fields**:
   - \( u \)-velocity: Strong horizontal flow near the lid, decreasing toward the bottom.
   - \( v \)-velocity: Vertical flow components forming a circulatory pattern.
2. **Pressure Field**: Smooth pressure distribution with gradients driving the flow.
3. **Streamlines**: A primary vortex in the cavity center, with possible secondary vortices in the corners for \( \text{Re} = 100 \).
4. **Animated GIF**: Shows the evolution of the \( u \)-velocity field over iterations, illustrating convergence.
5. **Convergence**: The continuity residual should drop below \( 10^{-5} \), indicating a physically consistent solution.

The results should match qualitative features of the lid-driven cavity flow, such as a central vortex and boundary layer effects near the walls.

---

## Suggestions for Improving the Project
1. **Grid Refinement**:
   - Increase the grid resolution (e.g., \( 100 \times 100 \)) to capture finer flow details, especially secondary vortices.
   - Implement adaptive mesh refinement near the walls for better accuracy.
2. **Higher Reynolds Numbers**:
   - Test higher \( \text{Re} \) (e.g., 400, 1000) to study transitional flows, adjusting the relaxation factors and tolerance as needed.
3. **Advanced Convection Schemes**:
   - Replace the upwind scheme with higher-order schemes (e.g., central differencing) for improved accuracy, balancing stability.
4. **Parallelization**:
   - Use FiPy's parallel capabilities or integrate with libraries like PETSc to speed up computations for larger grids.
5. **Additional Visualizations**:
   - Add vorticity contour plots to highlight rotational flow structures.
   - Include velocity vector plots to complement streamlines.
6. **Quantitative Validation**:
   - Compare results with benchmark data (e.g., Ghia et al., 1982) for \( u \)- and \( v \)-velocity profiles along the cavity centerlines.
7. **Transient Simulation**:
   - Modify the code to simulate transient flow by reducing the pseudo-timestep and tracking time evolution.
8. **Error Handling**:
   - Add checks for numerical stability (e.g., CFL condition) and graceful handling of non-convergence.
9. **Interactive Interface**:
   - Integrate with Jupyter Notebook or a GUI framework (e.g., Tkinter) to allow parameter adjustments and real-time visualization.

---

## Installation and Usage
1. **Clone the Repository**:
   ```bash
   git clone <repository_url>
   cd lid-driven-cavity-simulation
   ```
2. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
   Create a `requirements.txt` with:
   ```
   fipy==3.4.2
   numpy==1.26.4
   matplotlib==3.8.3
   imageio==2.34.0
   ```
3. **Run the Simulation**:
   ```bash
   python lid_driven_cavity.py
   ```
4. **Outputs**:
   - Contour plots: `u_velocity_field.png`, `v_velocity_field.png`, `pressure_field.png`.
   - Streamline plot: `stream_function_and_lines.png`.
   - Animated GIF: `lid_driven_cavity_evolution.gif` (if enabled).

---

## Contributing
Contributions are welcome! Please follow these steps:
1. Fork the repository.
2. Create a new branch (`git checkout -b feature-branch`).
3. Make changes and commit (`git commit -m "Add feature"`).
4. Push to the branch (`git push origin feature-branch`).
5. Open a pull request.

---

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

---

## Acknowledgments
- **FiPy Developers**: For providing a robust finite volume PDE solver.
- **CFD Community**: For extensive literature and benchmark data on the lid-driven cavity problem.
- **Ghia et al. (1982)**: For providing reference solutions for validation.

---

## Contact
For questions or suggestions, please open an issue on GitHub .

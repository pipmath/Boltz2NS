# Lattice BGK and Navier-Stokes Numerical Solvers (Boltz2NS)
MATLAB code to solve the Lattice BGK and Navier-Stokes equation, and compare the numerical solutions in order to determine the rate of the hydrodynamic limit of the macroscopic vorticity in 2D.

# Files
- `boltz_init.m`: Script initializes variables and workspace for running Navier-Stokes and LBGK Boltzmann solvers.
- `BoltzSolver_RK4.m`: Script used to solve the 2D LBGK Boltzmann, using a vorticity initial condition
- `dnsNS_IMEX.m`: Script used to solve the 2D Navier-Stokes, using the vorticity transport equation
- `Genplts.m`: Script runs Navier-Stokes solver (`dnsNS_IMEX.m`) and LBGK Boltzmann solver (`BoltzSolver_RK4.m`) using parameters given in `boltz_init.m`, and plots solutions for comparison 

# How to Use
Set parameters (viscosity, Knudesen number, spatial domain, time interval, spatial resolution, time step size, and number of time steps to save) in `boltz_init.m`, which is then used in the other scripts. `BoltzSolver_RK4.m` and `dnsNS_IMEX.m` can be ran individually and used to save files. Running `Genplts.m` ensures the necessary data exists for given parameters (and if not it runs the necessary solver), then plots the vorticity fields for four snapshots in time as well as a comparison of energy quantites between the Navier-Stokes and LBGK Boltzmann solutions.

The initial conditions *TG* and *TGvort* are already included in `boltz_init.m`, which represent the Taylor-Green and perturbed Taylor-Green, respectively. Since the LBGK problem is "stiff", the time step must be refined as the Knudsen number is decreased. 


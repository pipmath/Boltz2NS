# Lattice BGK and Navier-Stokes Numerical Solvers (Boltz2NS)
MATLAB code to solve the Lattice BGK Boltzmann (D2Q9) and Navier-Stokes equation, and compare the numerical solutions in order to determine the rate of the hydrodynamic limit of the macroscopic vorticity in 2D.

# Files
- `boltz_init.m`: Script initializes variables and workspace for running Navier-Stokes and LBGK Boltzmann solvers.
- `BoltzSolver_RK4.m`: Script used to solve the 2D LBGK Boltzmann, using a vorticity initial condition
- `dnsNS_IMEX.m`: Script used to solve the 2D Navier-Stokes, using the vorticity transport equation
- `Genplts.m`: Script runs Navier-Stokes solver (`dnsNS_IMEX.m`) and LBGK Boltzmann solver (`BoltzSolver_RK4.m`) using parameters given in `boltz_init.m`, and plots solutions for comparison 

# How to Use
Set parameters (viscosity, Knudesen number, spatial domain, time interval, spatial resolution, time step size, and number of time steps to save) in `boltz_init.m`, which is then used in the other scripts. `BoltzSolver_RK4.m` and `dnsNS_IMEX.m` can be ran individually and used to save files. Running `Genplts.m` ensures the necessary data exists for given parameters (and if not it runs the necessary solver), then plots the vorticity fields for four snapshots in time as well as a comparison of energy quantites between the Navier-Stokes and LBGK Boltzmann solutions.

The initial conditions *TG* and *TGvort* are already included in `boltz_init.m`, which represent the Taylor-Green and perturbed Taylor-Green, respectively. Since the LBGK problem is "stiff", the time step must be refined as the Knudsen number is decreased. 

# Citing
Work has been published in Nonlinearity. The paper can be found [here](https://doi.org/10.1088/1361-6544/adca81).

Zhongyang Gu, Xin Hu, Pritpal Matharu, Bartosz Protas, Makiko Sasada, and Tsuyoshi Yoneda. *The incompressible Navier–Stokes limit from the discrete-velocity BGK Boltzmann equation.* Nonlinearity **38**, 055014, 2025, https://doi.org/10.1088/1361-6544/adca81 

Bibtex:
```
@article{Nonlinearity38_2025,
doi = {10.1088/1361-6544/adca81},
url = {https://dx.doi.org/10.1088/1361-6544/adca81},
year = {2025},
month = {apr},
publisher = {IOP Publishing},
volume = {38},
number = {5},
pages = {055014},
author = {Gu, Zhongyang and Hu, Xin and Matharu, Pritpal and Protas, Bartosz and Sasada, Makiko and Yoneda, Tsuyoshi},
title = {{The incompressible Navier–Stokes limit from the discrete-velocity BGK Boltzmann equation}},
journal = {Nonlinearity},
abstract = {In this paper, we extend the Bardos–Golse–Levermore program (Bardos et al 1993 Commun. Pure Appl. Math. 46 667–753) to prove that a local weak solution to the d-dimensional incompressible Navier–Stokes equations () can be constructed by taking the hydrodynamic limit of a discrete-velocity Boltzmann equation with a simplified Bhatnagar–Gross–Krook collision operator. Moreover, in the case when the dimension is , we characterise the combinations of finitely many particle velocities and probabilities that lead to the incompressible Navier–Stokes equations in the hydrodynamic limit. Numerical computations conducted in two-dimensional indicate that in the case of the simplest velocity lattice (D2Q9), the rate with which this hydrodynamic limit is achieved is of order , where  is the Knudsen number. For the future investigations, it is worth considering if the hydrodynamic limit of the discrete-velocity Boltzmann equation can be also rigorously justified in the presence of non-trivial boundary conditions.}
}
```

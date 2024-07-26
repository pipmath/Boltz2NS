% BoltzSolver_RK4.m: Script used to solve the 2D LBGK Boltzmann, using a vorticity
% initial condition.
%
% Computes a vorticity using the LBGK with D2Q9 stencil for the (independent)
% velocity field, and a pseudo-spectral Galerkin approach with dealiasing using a
% Gaussian spectral filter along with the standard RK4 time stepping method.
% Computations are time stepped in Fourier space.
%
% Author: Pritpal 'Pip' Matharu
% Numerical Analysis, Department of Mathematics
% KTH Royal Institute of Technology
% Date: 2024/07/24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BoltzSolver_RK4
%% Initialize Solver
% Initialize parameters, solver matrices, and save arrays for LBGK Boltzmann 
[g_hat, wgts, saveBoltz, tvals, derivs, disct, params, T] = solve_init;
% Check if case has already been run and file is stored
if exist(params.str, 'file')
    prompt = " Data already exists. Do you want to continue and overwrite? [Y]: ";
    resp = input(prompt,"s");
    if strcmpi(resp, 'n')
        disp('Exiting Boltzmann solver!')
        return
    end
end

%% Main time stepping loop
for ti = 2:tvals.Nd+1

    % Assign data and use RK4 time stepping
    g_hat = RK4stepping(g_hat);
    % Update time
    T = T + disct.dt;

    % Saving the vorticity field and information
    if ( mod(ti, tvals.TSCREEN) == 1 )
        % Position for saving
        ipos = floor(ti/tvals.TSCREEN)+1;
        % Compute velocity components
        [u1_hat,  u2_hat, ~, rho] = rhou1u2(g_hat, derivs.Kx, derivs.Ky, derivs.Pois);
        % Compute vorticity field
        w  = real( ifft2( derivs.Kx.*u2_hat - derivs.Ky.*u1_hat ) );
        % Save current iteration of time
        saveBoltz.t(ipos)    = T;
        % Determine the enstrophy in the system and save
        saveBoltz.Enst(ipos) = 0.5 * sum( sum( w.*w ) ) * disct.dx * disct.dy;
        % Mass in system
        saveBoltz.M(ipos)          = sum(sum(rho));
        % Save vorticity for post-processing
        saveBoltz.w(:, :, ipos)  = w;
        % Save corresponding time value
        saveBoltz.tplt(ipos)= T;

        % Ensure the solution has not produced NaN, if so exit
        if ( max( max( isnan(w)) ) )
            disp(['Error with solution at t = ' num2str(T) ', NaN value']);
            return
        else
            fprintf(' t = %d, OK.  Enstrophy = %f \n', T, saveBoltz.Enst(ipos));
        end
    end

end
% Save time elapsed time
time_Boltz = toc;

% We want to save the final time
wBoltz     = w;

% Save Data
save(params.str, 'saveBoltz', 'wBoltz', 'time_Boltz')

disp('============================================================================');
fprintf(' Boltzmann solver complete. Clock Time: %f \n', time_Boltz);
disp('============================================================================');

%% ****************************** NESTED FUNCTIONS **********************************
    % RK4 Time Stepping
    function g_hat = RK4stepping(g_hat)

        % RK4 time stepping method
        k1 = disct.dt*RHSgeqn( g_hat );            % First step
        k2 = disct.dt*RHSgeqn( g_hat + 0.5*k1 );   % Second step
        k3 = disct.dt*RHSgeqn( g_hat + 0.5*k2 );   % Third step
        k4 = disct.dt*RHSgeqn( g_hat +     k3 );   % Fourth step
        g_hat = g_hat + (k1 + 2*k2 + 2*k3 + k4)/6; % Solution at T + dt

        % Solve the RHS of the Boltzmann equation
        function RHS = RHSgeqn(g_hat)
            % Compute density via summation, and velocity via weighted sums (indep v)
            [u1_hat,  u2_hat, rhoh] = rhou1u2( g_hat, derivs.Kx, derivs.Ky, ... 
                                               derivs.Pois );

            % Compute geq from velocity and rho
            geqh(:, :, :) = u2geq_hat( u1_hat, u2_hat, rhoh, wgts, params.ep, ...
                                       params.cs, derivs.filt );

            % RHS term of Boltzmann equation
            RHS(:,:,1:9) = derivs.FTh(:, :, 1:9).*g_hat(:,:,1:9) + geqh(:,:,1:9)...
                                                            /(params.ep^2*params.nu);
        end % End of RHSgeqn
    end % End of RK4stepping

end % End of BoltzSolver_RK4

%% ****************************** FUNCTIONS *****************************************
% -----------------------------------------------------------------------------------
% FUNCTION: solve_init
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2024/07/25
%
% Function to initialize solvers, parameters, and matrices for D2Q9 LBGK Boltzmann
%
% OUTPUT
% g_hat ...... Particle velocity distribution (with nonlinear effects)
% wgts ....... Weights for D2Q9 Lattice Boltzmann
% tsave ...... Arrays for saving in time (M, Enst, t, w, tplt)
% tvals ...... Values for time loop (Nd, TSCREEN)
% derivs ..... Spatial derivative operators (Kx, Ky, Pois, filt)
% disct ...... Discretization values (dt, dx, dy)
% params ..... Physical parameters (nu, ep)
% T0 ......... Initial Time
%
% FORMAT
% [g_hat, wgts, tsave, tvals, derivs, disct, params, T0] = solve_init
%
% -----------------------------------------------------------------------------------
function [g_hat, wgts, tsave, tvals, derivs, disct, params, T0] = solve_init

% Initialize system parameters, variables, and workspace
boltz_init; 

disp( ' Boltzmann Solver');
disp('============================================================================');

% Numerical weights for Lattice Boltzmann
wgts(1:4) = 1.0/9.0;  % Major directions
wgts(5:8) = 1.0/36.0; % Minor directions
wgts(9)   = 4.0/9.0;  % Zeroth direction
wgts      = wgts(:);

%% Transform from NS initial condition to Boltzmann
% Take FT of NS initial condition and compute velocity field from vorticity
w_0hat  =  fft2(w_0);
u01_hat = -( Ky./Pois_hat ).*w_0hat; % Velocity in x
u02_hat =  ( Kx./Pois_hat ).*w_0hat; % Velocity in y
% Arbitrary constant initial density (rho= 1)
rho0      = zeros(Nx,Ny);  % Preallocate
rho0(:,:) = 1/numel(rho0); % Normalization factor for FFT
rhoh      = fft2(rho0);    % FFT

% Compute geq (used as initial condition for Boltzmann) from velocity and rho
g_hat = u2geq_hat(u01_hat, u02_hat, rhoh, wgts, ep, c_s, dealias);

% Compute density via summation, and velocity via weighted sums (independent v)
[u1_hat,  u2_hat, ~, rho] = rhou1u2(g_hat, Kx, Ky, Pois_hat);
% Compute Mass (after LBGK approx) using Gaussian quadrature
M0  = sum(sum(rho));      % Mass in system

% Compute vorticity from velocity field
w_B = real( ifft2( u2_hat.*Kx - u1_hat.*Ky ) );
% Compute Enstrophy using Gaussian quadrature
E0  = 0.5 * sum( sum( w_B.*w_B ) ) * dx * dy; % Enstrophy of IC
fprintf(' Enstrophy of IC         = %f \n', E0);

% Compute advection term operator
FTh(:, :, 1) = (Kx);
FTh(:, :, 2) = (Ky);
FTh(:, :, 3) = (-Kx);
FTh(:, :, 4) = (-Ky);
FTh(:, :, 5) = (Kx+Ky);
FTh(:, :, 6) = (-Kx+Ky);
FTh(:, :, 7) = (-Kx-Ky);
FTh(:, :, 8) = (Kx-Ky);
FTh(:, :, 9) = 0;
derivs.FTh   = -( 1/(ep^2*nuB) + FTh/ep);

% Store spatial derivative operators
derivs.Kx   = Kx;       % matrix of wavenumbers in x direction
derivs.Ky   = Ky;       % matrix of wavenumbers in y direction
derivs.Pois = Pois_hat; % Operator to solve Poisson' equation
derivs.filt = dealias;  % Gaussian spectral filter for dealiasing

% Preallocate vectors
tsave.M    = zeros(Nsave+1,1);
tsave.Enst = zeros(Nsave+1,1);
tsave.t    = zeros(Nsave+1, 1);
tsave.w    = zeros(Nx, Ny, Nsave+1);
tsave.tplt = zeros(Nsave+1, 1);

% Store values at t=0
tsave.M(1)       = M0;  % Mass
tsave.Enst(1)    = E0;  % Enstrophy
tsave.t(1)       = T0;  % Time
tsave.w(:, :, 1) = w_B; % Vorticity Field
tsave.tplt(1)    = T0;  % Time for plots

% Store discretization values
disct.dt = dt; % Time step size
disct.dx = dx; % Spatial step size in x
disct.dy = dy; % Spatial step size in y

% Store values for time loop
tvals.Nd      = Nd;      % Number of time steps
tvals.TSCREEN = TSCREEN; % Screen update interval time

% Store physical parameters
params.nu  = nuB; % Viscosity for Boltzmann equation
params.ep  = ep;  % Knudsen number
params.cs  = c_s; % Speed of sound constant
% String for saving
params.str = sprintf('Boltz_%s_N%d_dt%de%d_T%de%d_NU%de%d_EP%de%d_Ns%d.mat', ...
    IC, N_DNS, dtfac, dtcoeff, Tfac, Tcoeff, nufac, nucoeff, epfac, epcoeff, Nsave);
end % End of solve_init

% -----------------------------------------------------------------------------------
% FUNCTION: u2geq_hat
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2024/07/25
%
% Function uses the velocity components, density, lattice weights, and epsilon to 
% determine the particle velocity distribution geq with nonlinear effects O(epsilon)
%
% INPUT
% u1_hat ...... Fluid velocity, x component
% u2_hat ...... Fluid velocity, y component
% rhoh ........ Fourier transform of particle density
% w ........... Weights for D2Q9 Lattice Boltzmann
% ep .......... Knudsen number
% c_s ......... Speed of sound constant
% dealias ..... Gaussian spectral filter for dealiasing
%
% OUTPUT
% g_hat ..... Particle velocity distribution (with nonlinear effects)
%
% FORMAT
% geqh = u2geq_hat(u1h, u2h, rhoh, w, ep, dealias)
%
% -----------------------------------------------------------------------------------
function geqh = u2geq_hat(u1h, u2h, rhoh, w, ep, c_s, dealias)

% Coefficients including speed of sound constant
C2 = 1/c_s^2;
C4 = 1/(2*c_s^4);

% Fourier transform of velocity components for nonlinear components
u1 = real(ifft2(u1h));
u2 = real(ifft2(u2h));

% Precompute products
u11  = u1.*u1;
u22  = u2.*u2;
u12  = u1.*u2;
uabs = (C2/2)*( u22 + u11 );

% Dealias nonlinear terms
u11h  = ep*fft2(u11).*dealias;
u22h  = ep*fft2(u22).*dealias;
u12h  = ep*fft2(u12).*dealias;
uabsh = ep*fft2(uabs).*dealias;
rhoh  = rhoh.*dealias;

% Compute Boltzmann initial condition from NS velocity field
geqh(:,:,1) = ( rhoh + C2*u1h          -uabsh + C4*u11h                       )*w(1);
geqh(:,:,2) = ( rhoh + C2*u2h          -uabsh + C4*u22h                       )*w(2);
geqh(:,:,3) = ( rhoh - C2*u1h          -uabsh + C4*u11h                       )*w(3);
geqh(:,:,4) = ( rhoh - C2*u2h          -uabsh + C4*u22h                       )*w(4);
geqh(:,:,5) = ( rhoh + C2*( u1h + u2h) -uabsh + C4*2*u12h + C4*u11h + C4*u22h )*w(5);
geqh(:,:,6) = ( rhoh + C2*(-u1h + u2h) -uabsh - C4*2*u12h + C4*u11h + C4*u22h )*w(6);
geqh(:,:,7) = ( rhoh + C2*(-u1h - u2h) -uabsh + C4*2*u12h + C4*u11h + C4*u22h )*w(7);
geqh(:,:,8) = ( rhoh + C2*( u1h - u2h) -uabsh - C4*2*u12h + C4*u11h + C4*u22h )*w(8);
geqh(:,:,9) = ( rhoh                   -uabsh                                 )*w(9);
end % End of u2geq_hat

% -----------------------------------------------------------------------------------
% FUNCTION: rhou1u2
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2024/07/25
%
% Compute density via summation, and velocity via weighted sums (independent var v)
%
% INPUT
% g_hat ........ Particle velocity distribution
% Kx ........... Matrix of wavenumbers for spatial derivative in x
% Ky ........... Matrix of wavenumbers for spatial derivative in y
% Pois_hat ..... Matrix operator to solve Poisson' equation
%
% OUTPUT
% u1_hat ..... Fluid velocity, x component
% u2_hat ..... Fluid velocity, y component
% rhoh ....... Fourier transform of particle density
% rho ........ Particle density
%
% FORMAT
% [u1_hat,  u2_hat, rhoh, rho] = rhou1u2(g_hat, Kx, Ky, Pois_hat)
%
% -----------------------------------------------------------------------------------
function [u1_hat,  u2_hat, rhoh, rho] = rhou1u2(g_hat, Kx, Ky, Pois_hat)
% Storage matrix
g        = zeros(size(g_hat));
% Compute inverse Fourier Transform of Boltzmann solution
g(:, :, 1:9) = real(ifft2(g_hat(:, :, 1:9)));

% Compute density via summation
rho  = sum(g, 3);
rhoh = fft2(rho); % FFT of density
% NS velocity components, computed via weighted sums in the independent variable v
u1_hat    = g_hat(:,:,1) - g_hat(:,:,3) + g_hat(:,:,5) ...
          - g_hat(:,:,6) - g_hat(:,:,7) + g_hat(:,:,8);
u2_hat    = g_hat(:,:,2) - g_hat(:,:,4) + g_hat(:,:,5) ...
          + g_hat(:,:,6) - g_hat(:,:,7) - g_hat(:,:,8);

% Leray projection for the velocity field
laP_hat = u1_hat.*Kx./Pois_hat + u2_hat.*Ky./Pois_hat;

u1_hat  = u1_hat - Kx.*laP_hat; % Velocity component in x
u2_hat  = u2_hat - Ky.*laP_hat; % Velocity component in x
end % End of rhou1u2

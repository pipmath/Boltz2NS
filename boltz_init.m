% boltz_init.m: Script initializes variables and workspace for running Navier-Stokes 
% and LBGK Boltzmann solvers.
%
% Author: Pritpal 'Pip' Matharu
% Numerical Analysis, Department of Mathematics
% KTH Royal Institute of Technology
% Date: 2024/07/24
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize Workspace
clear
close all
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% Switch parameters
% TG:     initial condition of the Taylor-Green vortices
% TGvort: initial condition of TG, with additional vortex in the middle
IC  = 'TGvort';

%% Initialize parameters for Boltzmann solver
% Adjustable physical parameters
nuB   = 1e-4;            % Corresponding viscosity for LBGK Boltzmann equation
ep    = 1.0;             % Knudsen number in LBGK Boltzmann equation
Lx    = 1.0;   Ly = Lx;  % Spatial domain x and y, using periodic box
T1    = 1.0;  T0 = 0.0;  % Final and initial Time
% Adjustable Numerical parameters
N_DNS = 128;             % Spatial resolution, in each direction
dt    = 2e-4;            % Time step size
Nsave = 100;             % Number of time steps saved

% Fixed physical parameters
nuNS    = nuB/3.0;              % Corresponding viscosity for Navier-Stokes
c_s     = (1/3)^(1/2);          % Speed of sound constant
% Fixed Numerical parameters
Nx      = N_DNS;  Ny = N_DNS;   % Number of discretization points in x and y
dx      = Lx/Nx;  dy = Ly/Ny;   % Spatial step size
Nd      = round(T1/dt);         % Number of time steps
TSCREEN = floor(T1/(dt*Nsave)); % Screen update interval time 

%% Setup wavenumbers
% Compute wavenumbers/ derivative operators in Fourier space
% x vector
KX_temp = (mod( (1:Nx)-ceil(Nx/2+1), Nx ) - floor(Nx/2));
KX_temp(round(Nx/2+1)) = round(Nx/2);       % Derivatives of even order
KX = (2*pi/Lx)*1i*ones(1, Ny)'*( KX_temp ); % matrix of wavenumbers in x direction
% y vector
KY_temp = (mod((1:Ny)' - ceil(Ny/2+1), Ny) - floor(Ny/2));
KY_temp(round(Ny/2+1)) = round(Ny/2);     % Derivatives of even order
KY = (2*pi/Lx)*1i*(KY_temp)*ones(1, Nx);  % matrix of wavenumbers in y direction
% Laplacian in Fourier space
Lap_hat        = KX.^2 + KY.^2;
% Operator to solve Poisson' equation
Pois_hat       = Lap_hat;
Pois_hat(1, 1) = 1; % Operator is inverted, so ensure the 0th wavenumber is constant

% For odd order of derivatives (to ensure the fft does not lose symmetry)
KX_temp(round(Nx/2+1)) = 0.0;               % Odd derivative orders (ensure symmetry)
KY_temp(round(Ny/2+1)) = 0.0;               % Odd derivative orders (ensure symmetry)
Kx = (2*pi/Lx)*1i*ones(1, Ny)'*( KX_temp ); % matrix of wavenumbers in x direction
Ky = (2*pi/Lx)*1i*(KY_temp)*ones(1, Nx);    % matrix of wavenumbers in y direction

% x, y grid for plotting
x  = linspace(0,Lx*(1-1/Nx),Nx)';
y  = linspace(0,Ly*(1-1/Ny),Ny)';
[Xg,Yg] = meshgrid(x,y);

% Determine the values of the fourier modes with cutoff
kmod    = sqrt( (abs(KX)).^2 + (abs(KY)).^2 );
kcut    = max(max(kmod))*(2/3);
% Use a Gaussian spectral filter for dealiasing
dealias = exp( -36 * (kmod/kcut).^36 );

%% Initial data
%% Determine the type of initial condition
% Obtain a initial condition for the vorticity (w0)
switch lower(IC)
    case {'tg'}
        % Taylor-Green Vortex
        a = 2;
        b = a;
        w_0 = 10*sin(2*pi*a*Xg).*sin(2*pi*b*Yg);
    case {'tgvort'}
        % Taylor-Green vortex with an additional vortex in the middle
        w_0 =(-sin(2*pi*Xg).*sin(2*pi*Yg))+ exp(-((Xg-0.5).^2+(Yg-0.5).^2)/(0.02));
    otherwise
        disp('Unknown Initial Condition!')
        return
end
%% Ensure initial condition has mean zero
% Determine Fourier transform
w_hat   = fft2( w_0 );
w_hat(1, 1) = 0.0; % Ensure zeroth Fourier mode is zero
% Transform back to physical space
w_0 = real(ifft2(w_hat));

%% For saving strings
% Epsilon value
epcoeff = floor(log10(ep));
epfac   = round(ep/(10^epcoeff));
% Viscosity Navier-Stokes value
nucoeff = floor(log10(nuNS));
nufac   = floor(nuNS/(10^nucoeff));
% Time step value
dtcoeff = floor(log10(dt));
dtfac   = floor(dt/(10^dtcoeff));
% Time value
Tcoeff = floor(log10(T1));
Tfac   = floor(T1/(10^Tcoeff));

% Display information
tic
disp('============================================================================');
disp(datetime);
disp([' Initial Condition: ' IC]);
disp( ' Parameters: ');
disp([' epsilon             = ' num2str(ep)]);
disp([' time step           = ' num2str(dt)]);
disp([' NS viscosity        = ' num2str(nuNS)]);
disp([' Spatial points (1D) = ' num2str(N_DNS)]);
disp([' Nsave               = ' num2str(Nsave)]);
disp([' T1                  = ' num2str(T1)]);
disp('============================================================================');


% dnsNS_IMEX.m: Script used to solve the 2D Navier-Stokes, using the  vorticity 
% transport equation
%
% Computes a DNS using a pseudo-spectral Galerkin approach with dealiasing
% using a Gaussian spectral filter and a globally third-order, four step
% IMEX time stepping method (which maintains good stability)
%
% Author: Pritpal 'Pip' Matharu
% Numerical Analysis, Department of Mathematics
% KTH Royal Institute of Technology
% Date: 2024/07/24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dnsNS_IMEX

%% Initialize Solver
% Initialize parameters
boltz_init; 
disp( ' Navier-Stokes Solver ');
disp('============================================================================');

% String for saving 
str  = sprintf('NS_%s_N%d_dt%de%d_T%de%d_NU%de%d_Ns%d.mat', ...
    IC, N_DNS, dtfac, dtcoeff, Tfac, Tcoeff, nufac, nucoeff, Nsave);
% Check if data exists
if exist(str, 'file')
    prompt = "Data already exists. Do you want to continue and overwrite? [Y]: ";
    resp = input(prompt,"s");
    if strcmpi(resp, 'n')
        disp('Exiting Navier-Stokes solver!')
        return
    end
end
% If continuing, save parameters and operator matrices
save(str, 'w_0', 'nuNS', 'Lx', 'Ly', 'Xg', 'Yg', 'Pois_hat', 'dealias', 'Kx', 'Ky', ...
    'N_DNS', 'dt', 'dx', 'dy', 'T1', 'str', '-v7.3')

%% Solve DNS of the 2D Navier-Stokes equation
dns_NS( w_0, Ky, Kx, Lap_hat, Pois_hat, dealias, nuNS, Ny, Nx, dt, dy, dx, ...
       TSCREEN, T0, T1, Nsave, str );
end



%% ****************************** FUNCTIONS *******************************
% -----------------------------------------------------------------------------------
% FUNCTION: u2geq_hat
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2024/07/25
%
% Function computes a DNS using a pseudo-spectral Galerkin approach with dealiasing
% using a Gaussian spectral filter and a globally third-order, four step IMEX time 
% stepping method (which maintains good stability)
%
% INPUT
% w0 ........... Initial vorticity
% Kx ........... Matrix of wavenumbers for spatial derivative in x
% Ky ........... Matrix of wavenumbers for spatial derivative in y
% Lap_hat ...... Laplacian in Fourier space 
% Pois_hat ..... Operator for Poisson's eqn 
% forc_hat ..... Forcing in Fourier space   
% dealias ...... Filtering matrix           
% nu ........... Viscosity for NS           
% Ny ........... Discretization points in y 
% Nx ........... Discretization points in x 
% dt ........... Step size in time          
% dy ........... Step size in y domain      
% dx ........... Step size in x domain      
% TSCREEN ...... Plotting update interval   
% T0 ........... Initial time               
% T1 ........... Final time                 
% Nsave ........ Number of entries to save  
% str .......... String for saving          
%
% OUTPUT
% w_NS ........ Vorticity, various times 
% Enst_NS ..... Enstrophy, all time points 
% t_NS ........ Time vector                
%
% FORMAT
% [w_NS, Enst_NS, t_NS] = dns_NS(w0, Ky, Kx, Lap_hat, Pois_hat, dealias, ...
%   nu, Ny, Nx, dt, dy, dx, TSCREEN, T0, T1, Nsave, str)
%
% -----------------------------------------------------------------------------------
function [w_NS, Enst_NS, t_NS] = dns_NS(w0, Ky, Kx, Lap_hat, Pois_hat, dealias, ...
           nu, Ny, Nx, dt, dy, dx, TSCREEN, T0, T1, Nsave, str)
% -----------------------------------------------------------------------------------
% Modify the Laplacian to include viscosity 
Lap_hat = nu*Lap_hat;

% Values for the IMEX method
% Coefficients from https://doi.org/10.1007/s10898-019-00855-1
% Journal of Global Optimization (2020) - Alimo, Cavaglieri, Beyhaghi, Bewley 
alphacoeff = dt *[343038331393.0/1130875731271.0 288176579239.0/1140253497719.0 ...
                  253330171251.0/677500478386.0  189462239225.0/1091147436423.0]';
betaIcoeff = dt *[35965327958.0/140127563663.0   19632212512.0/2700543775099.0  ...
                 -173747147147.0/351772688865.0  91958533623.0/727726057489.0]';
betaEcoeff = dt *[14.0/25.0                      777974228744.0/1346157007247.0 ...
                  251277807242.0/1103637129625.0 113091689455.0/220187950967.0]';
gammacoeff = dt *[0.0                           -251352885992.0/790610919619.0 ...
                -383714262797.0/1103637129625.0 -403360439203.0/1888264787188.0]';

% Empty storage arrays
t_NS       = zeros(Nsave+1, 1);
Enst_NS    = zeros(Nsave+1, 1);
w_NS       = zeros(Ny, Nx, Nsave+1);
tplt_NS    = zeros(Nsave+1, 1);

% Save the initial vorticity
w = w0;
w_NS(:, :, 1) = w;
% Transform vorticity to Fourier space
w_hat         = fft2(w);
% Determine the streamfunction by using: -lap(p) = w
psi_hat       = -w_hat ./ Pois_hat;

% Determine and save the enstrophy
% Compute L2 inner product using Gaussian quadrature
Enst_NS(1) = 0.5 * sum( sum( w.*w ) ) * dx * dy;

% Save initial time
t_NS(1)    = T0;
T          = T0;
tplt_NS(1) = T0;

% Display the mean of the initial condition (should be zero) and Enstrophy
fprintf('   Initial Enstrophy: %f \n   Mean: %f \n', Enst_NS(1), sum( sum( w ) ) * dx * dy);

% Matrix to be used in the IMEX method
conv0_hat(1:Ny, 1:Nx) = 0.0;

%% Main loop
for k = 2:round(T1/dt)+1
    % IMEX time stepping
    for rk = 1:4
        % Use computed stream function, and get the velocity and gradient of vorticity
        conv_hat = Jac_hat(psi_hat, w_hat, Ky, Kx); % Compute the Jacobian
        % Dealias and obtain correct sign for the Jacobian
        conv_hat = -conv_hat .* dealias; 

        % Compute solution at next sub-step using IMEX time-stepping method
        w2_hat = ( ...
                  ( 1.0 +  betaIcoeff(rk) * Lap_hat) .* w_hat + ...
                  betaEcoeff(rk) * conv_hat + ...
                  gammacoeff(rk) * conv0_hat ...
                 )./ ...
                 ( 1.0 - alphacoeff(rk) * Lap_hat );

        %% Update vorticity for next step
        w_hat     = w2_hat;

        % Update explicit part, for next substep
        conv0_hat = conv_hat;
        
        % Determine new streamfunction by using: -lap(p) = w
        psi_hat   = -w_hat ./ Pois_hat;
    end

    % Update time
    T = T + dt;
    
    %% Saving the vorticity field and information
    if ( mod(k, TSCREEN) == 1 )
        % Position for saving
        kpos = floor(k/TSCREEN)+1;
        % Save current iteration of time
        t_NS(kpos)    = T;
        % Transform back to real space
        w         = real(ifft2( w_hat ));
        % Determine the enstrophy in the system and save
        Enst_NS(kpos) = 0.5 * sum( sum( w.*w ) ) * dx * dy;

        % Save vorticity for post-processing
        w_NS(:, :, kpos)  = w;
        
        % Save corresponding time value
        tplt_NS(kpos)= T;
        
        % Ensure the solution has not produced NaN, if so exit
        if ( max( max( isnan(w)) ) )
            disp(['Error with solution at t = ' num2str(T) ', NaN value']);
            return
        else
            fprintf(' t = %d, OK.  Enstrophy = %f \n', T, Enst_NS(kpos));
        end
    end
    
    % Setting up the RHS vectors for the first substep
    conv0_hat(1:Ny, 1:Nx) = 0.0;    
end

fprintf('   Final Enstrophy: %f \n   Mean: %f \n', Enst_NS(kpos), sum( sum( w ) ) * dx * dy);

% Save time elapsed time
time_NS = toc;
% We want to save the final time as the new initial condition (starts at a statistically stationary point)
wNS     = w;
% Saving data
save(str, 'w_NS', 'wNS', 'Enst_NS', 't_NS', 'tplt_NS', 'time_NS', '-append')


disp('============================================================================');
fprintf(' Navier-Stokes solver complete. Clock Time: %f \n', time_NS);
disp('============================================================================');


end % End of dns_NS

% -----------------------------------------------------------------------------------
% FUNCTION: Jac_hat
%
% AUTHOR ... Pritpal Matharu
% DATE ..... 2023/12/07
%
% Compute the Jacobian determinant (convection term in 2D NS)
%               J = psi_x * w_y - psi_y * w_x 
% of f,g given in Fourier space.
%
% INPUT
% psi_hat ..... Streamfunction
% w_hat ....... Vorticity
% Ky .......... Matrix of wavenumbers for spatial derivative in y
% Kx .......... Matrix of wavenumbers for spatial derivative in x
%
% OUTPUT
% J ..... Jacobian determinant (convection term)
%
% FORMAT
% conv_hat = Jac_hat(psi_hat, w_hat, Ky, Kx)
%
% -----------------------------------------------------------------------------------
function conv_hat = Jac_hat(psi_hat, w_hat, Ky, Kx)
% Compute y derivative of stream function
u        = real(ifft2( Ky .* psi_hat ));
% Compute x derivative of stream function
v        = real(ifft2( Kx .* psi_hat ));

% Compute x derivative of vorticity
w_x      = real(ifft2( Kx .* w_hat ));
% Compute y derivative of vorticity
w_y      = real(ifft2( Ky .* w_hat ));
% Compute Jacobian, convective derivative (u,v).grad(w)
conv     = u.*w_x - v.*w_y;

% Transform to Fourier space
conv_hat = fft2( conv );

end % End of Jac_hat


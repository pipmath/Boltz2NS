% Genplts.m: cript used to compare solutions of the Navier-Stokes vs LBGK Boltzmann
%
% Author: Pritpal 'Pip' Matharu
% Numerical Analysis, Department of Mathematics
% KTH Royal Institute of Technology
% Date: 2024/07/24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize Solver
% Run Navier-Stokes Solver
dnsNS_IMEX; 
% Run LBGK Boltzmann Solver
BoltzSolver_RK4;

%% Plotting Navier-Stokes Vorticity
% Initialize parameters
boltz_init; 
% String for loading Navier-Stokes solver data 
strNS  = sprintf('NS_%s_N%d_dt%de%d_T%de%d_NU%de%d_Ns%d.mat', ...
                  IC, N_DNS, dtfac, dtcoeff, Tfac, Tcoeff, nufac, nucoeff, Nsave);
% Loading data
load(strNS)

% Positions for plotting comparisons
tref = length(tplt_NS);
tvals_plt = [1; ceil(tref/4); ceil(tref/2); tref];

% Plot Navier-Stokes Vorticity solutions
figure(1) 
clf
tfig = tiledlayout(2, 2);
for k = 1:length(tvals_plt)
    nexttile
    strT = sprintf('t=%d', tplt_NS(tvals_plt(k)));
    plt_contour(Xg, Yg, w_NS(:, :, tvals_plt(k)), strT)
end
strT = sprintf('Vorticity for %s, Navier-Stokes Solver', IC);
title(tfig, strT, 'Interpreter', 'latex')

%% Plotting Boltzmann Vorticity and Mass
% String for loading Navier-Stokes solver data 
strBoltz  = sprintf('Boltz_%s_N%d_dt%de%d_T%de%d_NU%de%d_EP%de%d_Ns%d.mat', IC, N_DNS, dtfac, dtcoeff, Tfac, Tcoeff, nufac, nucoeff, epfac, epcoeff, Nsave);
% Loading data
load(strBoltz)

% Positions for plotting comparisons
tref = length(saveBoltz.tplt);
tvals_plt = [1; ceil(tref/4); ceil(tref/2); tref];

% Plot Vorticity solutions
figure(2) 
clf
tfig = tiledlayout(2, 2);
for k = 1:length(tvals_plt)
    nexttile
    strT = sprintf('t=%d', saveBoltz.tplt(tvals_plt(k)));
    plt_contour(Xg, Yg, saveBoltz.w(:, :, tvals_plt(k)), strT)
end
strT = sprintf('Vorticity for %s, Boltzmann Solver', IC);
title(tfig, strT, 'Interpreter', 'latex')

% Mass
figure(3) 
clf
axfig = gca;
hold on
box on
set(axfig, 'yscale', 'log')
plt_t(saveBoltz.t, saveBoltz.M, 'Mass', 'Conservation of Mass, Boltzmann Solver', '--')

%% Compare Energy Quantities
figure(4) 
clf
tfig = tiledlayout(2, 2);
nexttile % panel 1
axfig = gca;
hold on
plt_t(t_NS, sqrt(2*Enst_NS), '$\|\omega(t)\|$', '$L^2$ Norm of Vorticity', '.-')
plt_t(saveBoltz.t, sqrt(2*saveBoltz.Enst), '$\|\omega(t)\|$', '$L^2$ Norm of Vorticity', '--')
legent{1} = sprintf('Navier-Stokes Solver ($\\omega$)');
legent{2} = sprintf('Boltzmann Solver ($\\omega^{\\varepsilon}$)');
legend(legent, 'location', 'northeast', 'Interpreter', 'Latex');

% Energy difference
nexttile % panel 2
axfig = gca;
hold on
box on
plt_t(t_NS, abs(sqrt(2*Enst_NS) - sqrt(2*saveBoltz.Enst)), '$|\|\omega(t)\| - \|\omega^{\varepsilon}(t)\||$', 'Difference of $L^2$ Norms of Vorticity', '--')
set(axfig, 'yscale', 'log')

% L2 Norm of Difference of Vorticity
w_diff            = abs( saveBoltz.w - w_NS ); 
Ediff = sqrt(sum( sum( w_diff(:, :, 1:end).*w_diff(:, :, 1:end) ) ) * dx * dy); 
Ediff = Ediff(:);
EnstNS = sqrt(sum( sum( w_NS(:, :, 1:end).*w_NS(:, :, 1:end) ) ) * dx * dy); 
EnstNS = EnstNS(:);

nexttile % panel 3
axfig = gca;
hold on
box on
plt_t(tplt_NS, Ediff, '$\|\omega(t) - \|\omega^{\varepsilon}(t)\|\|$', '$L^2$ Norm of Difference', '--')

% Relative Error
nexttile % panel 4
axfig = gca;
hold on
box on
plt_t(tplt_NS, Ediff./EnstNS, '$\|\omega(t) - \|\omega_{\varepsilon}(t)\|\|/\|\omega(t)\|$', 'Relative Difference', '--')
% Title
title(tfig, 'Comparison of Boltzmann Solver to Navier-Stokes Solver', 'Interpreter', 'latex')

%% ****************************** FUNCTIONS *******************************
% Function for plotting the vorticity field
function plt_contour(Xg, Yg, w, strT)
vmax = max(max(abs(w)));
vcont = linspace(-vmax,vmax,30);
[~, h] = contourf(Xg, Yg, w, vcont);
set(h,'LineColor','none');
cb=colorbar;
% cb.Position = cb.Position + [0 0 0 0];
% shading flat;colormap('jet');
str_title = sprintf(strT);
title(str_title, 'Interpreter', 'Latex')
axis square;
ylim([0 Yg(end, 1)]); xlim([0 Xg(1, end)])
xlabel('$x$', 'Interpreter','latex')
ylabel('$y$', 'Interpreter','latex')
end

%% ------------------------------------------------------------------------
% Function for plots dependent on time
function plt_t(t, f, strY, strT, line_sty)
% Set default line style
if nargin < 5
    line_sty = '.';
end
plot(t, f, line_sty)
ylabel(strY, 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title(strT, 'Interpreter', 'Latex')
xlim([t(1), t(end)])
end


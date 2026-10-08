
% Mean-flow preprocessor for the Incompact3d-IOIM TBL inlet.
% Copyright (c) 2026, Anjia Ying and Lin Fu. MIT License; see LICENSE.
%
% Nx, Ny, Lx, Ly, delta_in and dt must match TBL_case/input.i3d.
% Writes profiles/U_X*.txt and profiles/V_X*.txt for iin = 4.

clc
clear

addpath('./src/');

%% Model constants
kf = 0.344;              % skin-friction correlation constants (documented only;
Cfit = 1.770;            % the values actually used live in f_Round1)
imode = 1;               % 1: input Re_theta; 0: input Re_tau

%% DNS setting
%%% must be consistent with the settings of Xcompact3d in
%%% /DNS_implementation/Incompact3d-IOIM/TBL_case/input.i3d

Uinf = 1;                % free-stream velocity (Constant as 1 in Xcompact3d for TBLs)

Nx = 2401;
Ny = 301;

Lx = 545.6;
Ly = 50;

delta_in  = 5.23;

dt = 0.005;

dx = Lx / (Nx - 1);
xM = (0:Nx-1).' * dx;

beta = 1.55;

yp = f_generate_grid(Ny,beta,Ly);

%% Inlet state: choose nu so the composite profile hits the DNS d99
Re_tgt = 604.6;
d99_tgt = delta_in;

[~, ~, y_aux, retau_rec, retheta_rec, redelta_rec] = f_Round1(Re_tgt, imode);
nu = Uinf * d99_tgt / redelta_rec;

theta_init = d99_tgt/redelta_rec*retheta_rec;
Re_theta_init = Uinf * theta_init / nu;
Re_delta_init = redelta_rec;
Re_tau_init = retau_rec;

%% Allocate and march
m_U = zeros(Nx, Ny);
m_upinf = zeros(Nx, 1);
CfM = zeros(Nx, 1);
thetaM = zeros(Nx, 1);
Re_thetaM = zeros(Nx, 1);
Re_deltaM = zeros(Nx, 1);
Re_tauM = zeros(Nx, 1);
m_delta99 = zeros(Nx, 1);
m_deltav = zeros(Nx, 1);

thetaM(1) = theta_init;
Re_thetaM(1) = Re_theta_init;
Re_deltaM(1) = Re_delta_init;
Re_tauM(1) = Re_tau_init;

imode = 1;

[ue, uinf_plus, cf, ~, deltav, delta99, ~, ~, ~] = ...
    f_velocity(Ny, Re_theta_init,nu, yp, imode);
m_U(1, :) = ue * uinf_plus;
m_upinf(1) = uinf_plus;
m_delta99(1) = delta99;
m_deltav(1) = deltav;
CfM(1) = cf;

for iii = 2:Nx
    % Forward Euler on the momentum integral, Cf lagged from station iii-1
    % and then recomputed at the midpoint Reynolds number.
    thetaM(iii) = thetaM(iii-1) + cf / 2 * dx;
    Re_thetaM(iii) = Uinf * thetaM(iii) / nu;

    Re_theta05 = 0.5 * (Re_thetaM(iii) + Re_thetaM(iii-1));

    [ue, uinf_plus, cf, ~, deltav, delta99, ~, redelta_rec, retau_rec] = ...
        f_velocity(Ny, Re_theta05,nu, yp, imode);

    m_U(iii, :) = ue * uinf_plus;
    m_upinf(iii) = uinf_plus;
    CfM(iii) = cf;
    Re_tauM(iii) = retau_rec;
    Re_deltaM(iii) = redelta_rec;
    m_delta99(iii) = delta99;
    m_deltav(iii) = deltav;

    if mod(iii, 50) == 0 || iii == Nx
        fprintf('Marching %d / %d\n', iii, Nx);
    end
end

d99M = Re_deltaM ./ Re_thetaM .* thetaM;

%% Wall-normal velocity from continuity
dUdx = zeros(Nx, Ny);
for j = 1:Ny
    dUdx(:, j) = first_derivative(m_U(:, j) ./ m_upinf, xM);
end

m_V = zeros(Nx, Ny);
for iii = 1:Nx
    vm = 0;
    for j = 2:Ny
        vm = vm - 0.5 * (dUdx(iii, j) + dUdx(iii, j-1)) * (yp(j) - yp(j-1));
        m_V(iii, j) = vm;
    end
end
m_V = m_V .* m_upinf;

%% Calculate the parameters

% n_sp: Time interval for sampling (corresponding to dt_sp^+ = (n_sp*dt)^+ = 1)
% N_I: Samplings inside an IOIM block (corresponding to a time duration of N_I*n_sp*dt*U_infty/delta_{e,re} >= 10)
% N_Gamma: Length of the recycling pool. N_Gamma > N_I and (N_Gamma-N_I/2)*n_sp*dt*U_infty/delta_{e,re} >= 10
% N_averaging: largest update interval of the IOIM transfer function.
% input.i3d names this coefficient N_averaging; N_averaging >= 10*N_I.

x_re = 40*delta_in;
[~,ind_xre] = min(abs(xM - x_re));
utau = Uinf./m_upinf;
t_tau = m_deltav./utau;
n_sp    = floor(t_tau(ind_xre)/dt);                              % (n_sp*dt)^+ <= 1
N_I_min = ceil(10*m_delta99(ind_xre)/(n_sp*dt*Uinf));            % N_I*n_sp*dt*Uinf/delta99 >= 10
disp('Recommended bounds for IOIM-S in input.i3d:')
fprintf('  recycling station x/delta_in = %.3f, x = %.6g\n', xM(ind_xre)/delta_in, xM(ind_xre))
fprintf('  n_sp:        no larger than %d\n', n_sp)
fprintf('  N_I:         at least %d\n', N_I_min)
fprintf('  N_Gamma:     greater than N_I, and at least N_I/2 + %d\n', N_I_min)
fprintf('  N_averaging: at least 10*N_I\n')

%% Output data
% Copy profiles/ to DNS implementation/Incompact3d-IOIM/TBL_case/profiles.
% Do not replace input.i3d. iin = 4 requires U_X1.txt.
dir_name = 'profiles';
mkdir(dir_name)

for iii = 1:Nx

    fid = fopen(['./',dir_name,'/U_X',num2str(floor(iii)),'.txt'], 'w');
    fprintf(fid, '%.15E\n', m_U(iii,:)'./m_upinf(iii)*Uinf);
    fclose(fid);

    fid = fopen(['./',dir_name,'/V_X',num2str(floor(iii)),'.txt'], 'w');
    fprintf(fid, '%.15E\n', m_V(iii,:)'./m_upinf(iii)*Uinf);
    fclose(fid);
end

%% Cf and Re_theta against the DNS
% 
% fig = figure(1);
% set(fig, 'Position', [100, 100, 1000, 280]);
% ha = tight_subplot(1, 2, [0.125, 0.08], [0.22, 0.12], [0.10, 0.08]);
% 
% axes(ha(1));
% plot(xM / d99_tgt, CfM, 'r--', 'LineWidth', 1.0);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
% ylim([0, max(CfM) * 1.25]);
% xlim([0, max(xM) / d99_tgt]);
% xlabel('$(x-x_0)/\delta_{e,x_0}$', 'Interpreter', 'latex');
% ylabel('$C_f$', 'Interpreter', 'latex');
% ax = gca;
% ax.YAxis.Exponent = -3;
% grid on
% subtitle('$(a)$', 'Interpreter', 'latex', ...
%     'Position', [-10, max(CfM) * 1.25, 0], 'FontSize', 15);
% 
% axes(ha(2));
% plot(xM / d99_tgt, Re_thetaM, 'r--', 'LineWidth', 1.0);
% set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
% ylim([0, max(Re_thetaM) * 1.25]);
% xlim([0, max(xM) / d99_tgt]);
% xlabel('$(x-x_0)/\delta_{e,x_0}$', 'Interpreter', 'latex');
% ylabel('$Re_\theta$', 'Interpreter', 'latex');
% ax = gca;
% ax.YAxis.Exponent = 2;
% grid on
% subtitle('$(b)$', 'Interpreter', 'latex', ...
%     'Position', [-10, max(Re_thetaM) * 1.25, 0], 'FontSize', 15);

%% Local functions
function idx = nearest_index(reference, target)
    [~, idx] = min(abs(reference(:) - target));
end

function e = rel_l2(prediction, reference, weight)
    % Weighted relative L2 norm. prediction and reference are columns.
    prediction = prediction(:);
    reference = reference(:);
    weight = weight(:);
    num = sqrt(sum((prediction - reference).^2 .* weight) / sum(weight));
    den = sum(reference .* weight) / sum(weight);
    e = num / den;
end

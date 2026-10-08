%PREDICT_TBL_MEAN_EVOLUTION
% Streamwise march of an incompressible ZPG turbulent boundary layer.
% Copyright (c) 2026, Anjia Ying and Lin Fu. MIT License; see LICENSE.
%
% The mean profile at each station is a composite law
% (van Driest inner layer + Coles-type wake) whose wake factor is chosen
% so that Cf matches
%
%     Cf = 2 / (log(Re_theta)/kf + C)^2,   kf = 0.344, C = 1.770.
%
% Stations are coupled by the von Karman momentum integral for zero
% pressure gradient,
%
%     d(theta)/dx = Cf/2,
%
% integrated with a forward Euler step and a midpoint evaluation of Cf.
% The wall-normal velocity follows from continuity,
%
%     V(y) = -int_0^y dU/dx dy'.
%
% Reference data are the UPM/Torroja ZPG boundary-layer DNS
% (Sillero, Jimenez & Moser, Phys. Fluids 25, 105102, 2013),
% https://torroja.dmt.upm.es/turbdata/blayers/high_re/
%
% Expected MAT files (see README):
%   data/tbl_xdata_structure.mat     variable "val"  -> streamwise records
%   data/structure_all_profiles.mat  variable "val"  -> profile structures
%
% Outputs are written to figures/ and are not required to re-run.

%% Paths
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'src'));

data_dir = fullfile(root_dir, 'data');
fig_dir = fullfile(root_dir, 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

xdata_file = fullfile(data_dir, 'tbl_xdata_structure.mat');
prof_file = fullfile(data_dir, 'structure_all_profiles.mat');
if ~isfile(xdata_file) || ~isfile(prof_file)
    error('predict_tbl_mean_evolution:MissingData', ...
        ['Place tbl_xdata_structure.mat and structure_all_profiles.mat in %s.\n' ...
         'See README.md for the required fields.'], data_dir);
end

%% User parameters
x_start_over_d99 = 10;   % march starts at the first grid point with x >= this * d99(1)
Uinf = 1;                % free-stream velocity used to form Reynolds numbers
kf = 0.344;              % skin-friction correlation constants (documented only;
Cfit = 1.770;            % the values actually used live in f_Round1)
imode = 1;               % 1: close the profile on Re_theta

fprintf('Skin-friction correlation constants: kf = %.3f, C = %.3f.\n', kf, Cfit);

%% Load the DNS reference
xdata = load(xdata_file);
allx = xdata.val;

pdata = load(prof_file);
profile_uvw = pdata.val;
n_points = numel(profile_uvw);

% Station index of each stored profile: mean of the nearest matches in
% d99, Re_tau, Re_theta, displacement thickness and u_tau.
ind_points = zeros(n_points, 1);
for iii = 1:n_points
    ind1 = nearest_index(allx.d99, profile_uvw(iii).d99);
    ind2 = nearest_index(allx.Retau, profile_uvw(iii).Retau);
    ind3 = nearest_index(allx.Reth, mean(profile_uvw(iii).Reth));
    ind4 = nearest_index(allx.dels, profile_uvw(iii).dels);
    ind5 = nearest_index(allx.utau, profile_uvw(iii).utau);
    ind_points(iii) = floor((ind1 + ind2 + ind3 + ind4 + ind5) / 5);
end

fig = figure(10);
plot(allx.xgrid(1:end-2), allx.cf, 'b-');
hold on
plot(allx.xgrid(ind_points), allx.cf(ind_points), 'p');
xlabel('x');
ylabel('C_f');
title('DNS skin friction and matched profile stations');
grid on
exportgraphics(fig, fullfile(fig_dir, 'stations_cf.pdf'), 'ContentType', 'vector');
savefig(fig, fullfile(fig_dir, 'stations_cf.fig'));

ind_start = find(allx.xgrid >= x_start_over_d99 * allx.d99(1), 1, 'first');

%% Grid of the march
Nx = (numel(allx.xgrid) - 2) - ind_start + 1;
Lx = allx.xgrid(end-2) - allx.xgrid(ind_start);
dx = Lx / (Nx - 1);
xM = (0:Nx-1).' * dx;

% Wall-normal grid taken from the first stored profile (outer units).
yp = profile_uvw(1).yd99(:) * profile_uvw(1).d99;
Ny = numel(yp);

%% Inlet state: choose nu so the composite profile hits the DNS d99
Re_tgt = allx.Reth(ind_start);
d99_tgt = allx.d99(ind_start);

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
    theta05 = 0.5 * (thetaM(iii) + thetaM(iii-1));

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

%% Colors used for the figures
m_color1 = [0.215686274509804,0.494117647058824,0.721568627450980];
m_color2 = [0.894117647058824,0.101960784313725,0.109803921568627];

%% Figure 1: von Karman residual of the DNS reference itself
xN = allx.xgrid(1:end-2) - allx.xgrid(ind_start);

dxve = allx.cf.' / 2 - first_derivative(allx.theta.', allx.xgrid(1:end-2).');
ind_end = find(xN / d99_tgt >= 85, 1, 'first');

fig = figure(3);
set(fig, 'Position', [100, 100, 600, 300]);
plot(xN / d99_tgt, dxve ./ (allx.cf / 2).', 'LineWidth', 1,'Color',m_color1);
hold on
xline(xN(1) / d99_tgt, 'k:', 'LineWidth', 1);
xline(xN(end) / d99_tgt, 'k:', 'LineWidth', 1);
xline(xN(ind_start) / d99_tgt, 'k--', 'LineWidth', 1);
if ~isempty(ind_end)
    xline(xN(ind_end) / d99_tgt, 'k--', 'LineWidth', 1);
end
yticks([-1, -0.5, 0, 0.5, 1]);
ylim([-1.2, 1.2]);
xlim([xN(1) - 20, xN(end) + 20] / d99_tgt);
grid on
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
xlabel('$(x-x_0)/\delta_{e,x_{e.0}}$', 'Interpreter', 'latex');
ylabel('$\sigma_\varepsilon$', 'Interpreter', 'latex');
exportgraphics(fig, fullfile(fig_dir, 'test_vk.pdf'), 'ContentType', 'vector');
savefig(fig, fullfile(fig_dir, 'test_vk.fig'));

fprintf('Figures written to %s\n', fig_dir);

%% Figure 2: Cf and Re_theta against the DNS
fig = figure(1);
set(fig, 'Position', [100, 100, 1000, 280]);
ha = tight_subplot(1, 2, [0.125, 0.08], [0.22, 0.12], [0.10, 0.08]);

axes(ha(1));
plot(xN / d99_tgt, allx.cf, '-', 'LineWidth', 1.2,'Color',m_color1);
hold on
plot(xM / d99_tgt, CfM, '--', 'LineWidth', 1.2,'Color',m_color2);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
legend('DNS', 'Predicted', 'Interpreter', 'latex', 'Location', 'southwest');
ylim([0, max(CfM) * 1.25]);
xlim([0, max(xN) / d99_tgt]);
xlabel('$(x-x_0)/\delta_{e,x_0}$', 'Interpreter', 'latex');
ylabel('$C_f$', 'Interpreter', 'latex');
ax = gca;
ax.YAxis.Exponent = -3;
grid on
subtitle('$(a)$', 'Interpreter', 'latex', ...
    'Position', [-10, max(CfM) * 1.25, 0], 'FontSize', 15);

axes(ha(2));
plot(xN / d99_tgt, allx.Reth, '-', 'LineWidth', 1.2,'Color',m_color1);
hold on
plot(xM / d99_tgt, Re_thetaM, '--', 'LineWidth', 1.2,'Color',m_color2);
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 15);
ylim([0, max(Re_thetaM) * 1.25]);
xlim([0, max(xN) / d99_tgt]);
xlabel('$(x-x_0)/\delta_{e,x_0}$', 'Interpreter', 'latex');
ylabel('$Re_\theta$', 'Interpreter', 'latex');
ax = gca;
ax.YAxis.Exponent = 2;
grid on
subtitle('$(b)$', 'Interpreter', 'latex', ...
    'Position', [-10, max(Re_thetaM) * 1.25, 0], 'FontSize', 15);

exportgraphics(fig, fullfile(fig_dir, 'PredictionM_Cf_Re.pdf'), 'ContentType', 'vector');
savefig(fig, fullfile(fig_dir, 'PredictionM_Cf_Re.fig'));

%% Figure 3: mean U and V profiles, staggered on a log axis
fig = figure(2);
set(fig, 'Position', [60, 80, 1200, 550]);
ha = tight_subplot(2, 6, [0.12, 0.05], [0.12, 0.12], [0.06, 0.012]);

dy = (yp(3:end) - yp(1:end-2)) / 2;
dy = [dy(1); dy; dy(end)];
err_U = zeros(n_points, 1);
err_V = zeros(n_points, 1);
tags = {'$(a1)$','$(a2)$','$(a3)$','$(a4)$','$(a5)$','$(a6)$', ...
        '$(b1)$','$(b2)$','$(b3)$','$(b4)$','$(b5)$','$(b6)$'};

ymaxU = 0;
ymaxV = 0;
for iii = 1:n_points
    iix = ind_points(iii) - ind_start + 1;
    if iix < 1 || iix > Nx
        continue
    end
    ind_99 = find(yp >= 0.99 * profile_uvw(iii).d99, 1, 'first');
    ymaxU = max(ymaxU, max(profile_uvw(iii).ua(1:ind_99)));
    ymaxV = max(ymaxV, max(profile_uvw(iii).va(1:ind_99)));
end

for iii = 1:n_points
    iix = ind_points(iii) - ind_start + 1;
    if iix < 1 || iix > Nx
        continue
    end
    ind_99 = find(yp >= 0.99 * profile_uvw(iii).d99, 1, 'first');
    x = yp(1:ind_99) / m_deltav(iix);

    axes(ha(iii));
    semilogx(x, profile_uvw(iii).ua(1:ind_99), '-', 'LineWidth', 1.2,'Color',m_color1); hold on
    semilogx(x, m_U(iix, 1:ind_99), '--', 'LineWidth', 1.2,'Color',m_color2);
    err_U(iii) = rel_l2(m_U(iix, 1:ind_99).', profile_uvw(iii).ua(1:ind_99), dy(1:ind_99));
    xlim([x(2)*0.5, x(end)*1.4]);
    xticks([1, 10, 100]);
    ylim([0, ymaxU*1.08]);
    grid on
    set(gca, 'XMinorGrid', 'off', 'TickLabelInterpreter', 'latex', 'FontSize', 16);
    if iii == 1
        ylabel('$U_1^+$', 'Interpreter', 'latex');
    else
%         set(gca, 'YTickLabel', []);
    end
%     set(gca, 'XTickLabel', []);
    text(-0.32,1.1, tags{iii}, 'Units', 'normalized', ...
        'Interpreter', 'latex', 'FontSize', 16);

    axes(ha(iii+6));
    semilogx(x, profile_uvw(iii).va(1:ind_99), '-', 'LineWidth', 1.2,'Color',m_color1); hold on
    semilogx(x, m_V(iix, 1:ind_99), '--', 'LineWidth', 1.2,'Color',m_color2);
    err_V(iii) = rel_l2(m_V(iix, 1:ind_99).', profile_uvw(iii).va(1:ind_99), dy(1:ind_99));
    xlim([x(2)*0.5, x(end)*1.4]);
    xticks([1, 10, 100]);
    ylim([0, ymaxV*1.15]);
    yticks([0, 0.02, 0.04]);
    grid on
    set(gca, 'XMinorGrid', 'off', 'TickLabelInterpreter', 'latex', 'FontSize', 16);
    xlabel('$y^+$', 'Interpreter', 'latex');
    if iii == 1
        ylabel('$U_2^+$', 'Interpreter', 'latex');
    else
%         set(gca, 'YTickLabel', []);
    end
    text(-0.32,1.1, tags{iii+6}, 'Units', 'normalized', ...
        'Interpreter', 'latex', 'FontSize', 16);
end

legend(ha(1), {'DNS','Predicted'}, 'Interpreter', 'latex', ...
    'Orientation', 'horizontal', 'Position', [0.4,0.95,0.2,0.05]);

exportgraphics(fig, fullfile(fig_dir, 'PredictionM_MeanUV.pdf'), 'ContentType', 'vector');
savefig(fig, fullfile(fig_dir, 'PredictionM_MeanUV.fig'));
fprintf('Relative L2 errors (U, V) at the %d profile stations:\n', n_points);
disp([err_U, err_V]);

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

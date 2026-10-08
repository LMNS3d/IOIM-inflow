function [ue, uinf_plus, y_aux, retau, retheta, redelta] = ...
    f_Round1(Re_tgt, imode)
%F_ROUND1 Wake factor on a self-contained auxiliary grid.
%   [UE, UINF_PLUS, Y_AUX, RETAU, RETHETA, REDELTA] = F_ROUND1(RE_TAU,
%   RE_THETA, IMODE) is the only routine that builds the auxiliary
%   wall-normal grid and the only routine that iterates the wake factor.
%
%   Inputs are only the two Reynolds numbers and the closure switch:
%       IMODE = 0  Re_tau prescribed, Re_theta returned
%       IMODE = 1  Re_theta prescribed, Re_tau returned
%   No physical grid and no target thickness are accepted. The auxiliary
%   grid is built at delta99 = 1, i.e. Y_AUX is y/delta99. The composite
%   profile depends on y/delta99 only, so the physical scale is applied
%   later by F_VELOCITY.
%
%   The wake factor Pi is a secant root of
%       Cf_corr/Cf_profile - 1 = 0,
%       Cf_profile = 2/(U_e/u_tau)^2,
%       Cf_corr    = 2/(log(Re_theta)/kf + C)^2,
%   with kf = 0.344 and C = 1.770. Each trial calls F_SCHEME_GV.
%
%   After Pi has converged, the grid is extended from 1 to 2 (in units of
%   delta99) with the same clipped slope used previously inside
%   F_SCHEME_GV, so F_VELOCITY can interpolate a mesh taller than delta99.
%   UE and Y_AUX are that extended auxiliary solution. UE is U/U_e.

ny = 10000;
d99_hat = 1;
y_core = (0:ny) * (d99_hat / ny);

PI_0 = 0.99;
PI_1 = 1.0;
cf_tol = 1e-10;

if imode == 0
    Re_tau = Re_tgt;
    Re_theta = 0;
elseif imode == 1
    Re_tau = 0;
    Re_theta = Re_tgt;
end

[~, uinf_plus, retau, retheta, redelta, err_0] = ...
    profile_residual(y_core, Re_tau, Re_theta, imode, PI_0);
[~, ~, ~, ~, ~, err_1] = ...
    profile_residual(y_core, Re_tau, Re_theta, imode, PI_1);

while true
    factor = 1.0;
    while true
        PI = PI_1 - factor * err_1 / ((err_1 - err_0) / (PI_1 - PI_0));
        if PI > PI_1 / 2
            break
        end
        factor = factor / 2;
    end

    [ue_core, uinf_plus, retau, retheta, redelta, err] = ...
        profile_residual(y_core, Re_tau, Re_theta, imode, PI);

    if abs(err) <= cf_tol || PI_1 <= 1e-5
        break
    end

    PI_0 = PI_1;
    PI_1 = PI;
    err_0 = err_1;
    err_1 = err;
end

[ue, y_aux] = extend_auxiliary(y_core, ue_core);
end

function [ue, uinf_plus, retau, retheta, redelta, err] = ...
    profile_residual(y_core, Re_tau, Re_theta, imode, PI_wake)
% One F_SCHEME_GV evaluation and the skin-friction closure residual.
kf = 0.344;
C = 1.770;

[ue, uinf_plus, retau, retheta, redelta] = ...
    f_Scheme_GV(y_core, Re_tau, Re_theta, imode, PI_wake);

Cf_profile = 2 / uinf_plus^2;
Cf_corr = 2 / (log(retheta) / kf + C)^2;
err = Cf_corr / Cf_profile - 1;
end

function [ue_ext, y_ext] = extend_auxiliary(y_core, ue_core)
% Extend from delta99 to 2*delta99. Slope is taken from the last interval
% of the core profile and reduced by 1.5, then clipped at the free stream.
ny = numel(y_core) - 1;
dy = y_core(2) - y_core(1);
y_ext = (0:(2 * ny)) * dy;

ue_ext = [ue_core, zeros(1, ny)];
dudy = (ue_ext(ny+1) - ue_ext(ny)) / (y_ext(ny+1) - y_ext(ny));
for j = ny+2:numel(y_ext)
    ue_ext(j) = ue_ext(ny+1) + dudy * (y_ext(j) - y_ext(ny+1)) / 1.5;
    ue_ext(j) = min(ue_ext(j), 1);
end
end

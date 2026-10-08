function [ue, uinf_plus, retau, retheta, redelta] = ...
    f_Scheme_GV(y_aux, retau, retheta, imode, PI_wake)
%F_SCHEME_GV Composite profile on a supplied auxiliary grid.
%   [UE, UINF_PLUS, RETAU, RETHETA, REDELTA] = F_SCHEME_GV(Y_AUX, RETAU,
%   RETHETA, IMODE, PI_WAKE) builds
%
%       U+(y+) = U_inner+(y+) + PI_WAKE * U_wake+(y+)
%
%   on the auxiliary wall-normal grid Y_AUX and iterates the unknown
%   Reynolds number until the momentum thickness of that profile is
%   consistent with the prescribed one.
%
%   Y_AUX must cover [0, delta99], with Y_AUX(1) = 0 and Y_AUX(end) =
%   delta99. This function does not build that grid, does not extend it
%   past delta99, and does not interpolate onto the physical mesh.
%   Those steps belong to F_ROUND1 and F_VELOCITY.
%
%   IMODE = 0  Re_tau is prescribed; Re_theta is returned.
%   IMODE = 1  Re_theta is prescribed; Re_tau is returned.
%
%   UE is U/U_e on Y_AUX, normalized so that U(y = delta99) = 0.99*U_e.
%   UINF_PLUS is U_e/u_tau. The profile shape in y/delta99 does not
%   depend on the value of delta99; only the product Re_tau * y/delta99
%   enters the wall law.

y_aux = y_aux(:).';
d99 = y_aux(end);
if y_aux(1) ~= 0 || d99 <= 0
    error('f_Scheme_GV:BadGrid', ...
        'Y_AUX must start at 0 and end at a positive delta99.');
end

if imode == 0
    retheta = 1000;   % initial guess; overwritten by the iteration
else
    retau = 200;
end

absdiff = inf;
niter = 0;
max_iter = 10000;
tol = 1e-5;
n = numel(y_aux);

while absdiff > tol
    niter = niter + 1;
    if niter > max_iter
        warning('f_Scheme_GV:MaxIter', ...
            'GFM-VD iteration stopped at %d steps (residual %.3e).', ...
            max_iter, absdiff);
        break
    end

    if imode == 0
        retheta_old = retheta;
    else
        retau_old = retau;
    end

    yp = (y_aux / d99) * retau;
    [up_i, up_o] = f_ITBL_mean(yp(:), retau);
    up_r = up_i + PI_wake * up_o;

    % Normalize so the 99% thickness is exactly the prescribed d99.
    uinf_plus = up_r(end) / 0.99;
    ue = (up_r / uinf_plus).';

    th = 0;
    for j = 2:n
        th = th + 0.5 * (ue(j) * (1 - ue(j)) + ue(j-1) * (1 - ue(j-1))) ...
            * (y_aux(j) - y_aux(j-1));
    end

    if imode == 0
        retheta = retau * th / d99 * uinf_plus;
        absdiff = abs(retheta - retheta_old);
    else
        retau = retheta * d99 / th / uinf_plus;
        absdiff = abs(retau - retau_old);
    end
end

redelta = retheta * d99 / th;
end

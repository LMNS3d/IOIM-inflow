function [UpN, UpL] = f_ITBL_mean(yp, Ret)
%F_ITBL_MEAN Incompressible composite mean profile on a wall-normal grid.
%   [UPN, UPL] = F_ITBL_MEAN(YP, RET) returns the inner (van Driest) and
%   wake contributions to the mean velocity, both in wall units, sampled
%   at the wall-normal stations YP (also in wall units). RET is the
%   friction Reynolds number used to scale the wake.
%
%   Inner contribution
%       d+ = (1 - exp(-y+/A+))^2,   A+ = 17
%       nu_t+ = kappa * y+ * d+,    kappa = 0.41
%       dU+/dy+ = 1 / (1 + nu_t+)
%   integrated from the wall with the trapezoidal rule.
%
%   Wake contribution (Coles-type, sine distribution of the wake shear)
%       dU_wake+/dy = (pi/kappa) / Re_tau * sin(pi * y/delta)
%   integrated from the wall. The caller multiplies UPL by a wake factor
%   Pi and adds it to UPN.
%
%   Stations above y+ = Re_tau are clipped to the boundary-layer edge
%   before interpolation, then the edge value is held constant.

Aplus = 17;
kappa = 0.41;

yp = yp(:);
ny = max(min(ceil(2 * Ret), 10000), numel(yp));
y = linspace(yp(1), Ret, ny)';

eta = min(y / Ret, 1);
ypl = eta * Ret;
damping = (1 - exp(-ypl / Aplus)).^2;
nut = kappa .* ypl .* damping;

% Inner law: integrate 1/(1+nu_t+) in y+.
U1 = zeros(ny, 1);
for j = 2:ny
    uint = 0.5 * (1 / (1 + nut(j)) + 1 / (1 + nut(j-1)));
    U1(j) = U1(j-1) + uint * (ypl(j) - ypl(j-1));
end

% Wake shear: (pi/kappa)/Re_tau * sin(pi*eta), integrated in y.
S = (1 / Ret) * (pi / kappa) * sin(pi * eta);
U2 = zeros(ny, 1);
for j = 2:ny
    U2(j) = U2(j-1) + 0.5 * (S(j) + S(j-1)) * (y(j) - y(j-1));
end

% Hold the edge value so interp1 never returns NaN above the mesh.
y = [y; yp(end) + 1e-10];
U1 = [U1; U1(end)];
U2 = [U2; U2(end)];

UpN = interp1(y, U1, yp, 'linear', 0);
UpL = interp1(y, U2, yp, 'linear', 0);
end

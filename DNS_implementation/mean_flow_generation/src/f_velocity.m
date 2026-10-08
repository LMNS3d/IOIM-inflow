function [ue, uinf_plus, cf, th, deltav, delta99, ...
          retheta_rec, redelta_rec, retau_rec] = ...
    f_velocity(ny, Re_tgt,nu, y, imode)
%F_VELOCITY Interpolate the auxiliary profile onto the physical grid.
%   [UE, UINF_PLUS, CF, TH, DELTAV, DELTA99, RETHETA, REDELTA, RETAU] =
%   F_VELOCITY(NY, THETA_TGT, RE_TAU, RE_THETA, Y, IMODE) maps the
%   auxiliary solution from F_ROUND1 onto the physical wall-normal grid Y.
%
%   This routine does not build the auxiliary grid and does not iterate
%   the wake factor. Both of those steps are inside F_ROUND1, which is
%   called once. The profile shape is a function of y/delta99 only, so
%   changing delta99 does not require a new wake-factor solve.
%
%   What remains here is the physical scale. DELTA99 is chosen so that the
%   momentum thickness integrated on Y equals THETA_TGT to a relative
%   tolerance of 1e-3. UE on Y is the linear interpolant of the extended
%   auxiliary profile; samples beyond the auxiliary mesh are held at 1,
%   matching the previous edge treatment.
%
%   IMODE is forwarded to F_ROUND1 (0: Re_tau given, 1: Re_theta given).
%   NY must equal numel(Y). The output positions are unchanged from v1:
%   the sixth output is the 99% thickness, not a momentum thickness.

y = y(:).';
if numel(y) ~= ny
    error('f_velocity:BadSize', 'NY must equal the length of Y.');
end

[ue_aux, uinf_plus, y_hat, retau_rec, retheta_rec, redelta_rec] = ...
    f_Round1(Re_tgt, imode);

U_e = 0.99;
delta99 = nu*redelta_rec/U_e;

ue = interpolate_to_physical(y, y_hat * delta99, ue_aux);

cf = 2 / uinf_plus^2;
deltav = delta99 / retau_rec;
th = retheta_rec * delta99 / redelta_rec;
end

function ue = interpolate_to_physical(y, y_aux, ue_aux)
% Linear map from the extended auxiliary mesh. Outside it, hold Ue = 1.
ue = ones(1, numel(y));
inside = y < y_aux(end);
ue(inside) = interp1(y_aux, ue_aux, y(inside), 'linear');
end

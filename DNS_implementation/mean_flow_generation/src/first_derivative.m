function dfdx = first_derivative(f, x)
%FIRST_DERIVATIVE First derivative of a 1-D sample by finite differences.
%   DFDX = FIRST_DERIVATIVE(F, X) returns dF/dX on the same nodes as F.
%   Interior nodes use a second-order central difference. The first node
%   uses a first-order forward difference and the last node a first-order
%   backward difference.
%
%   F and X must be vectors of equal length. The output is a column vector.
%
%   This is a non-uniform-grid scheme: the local spacing is taken from X,
%   so a stretched wall-normal mesh is admissible. It is not suitable for
%   periodic data.

f = f(:);
x = x(:);

if numel(f) ~= numel(x)
    error('first_derivative:SizeMismatch', ...
        'F and X must have the same number of elements.');
end
if numel(f) < 2
    error('first_derivative:TooShort', ...
        'At least two samples are required.');
end

n = numel(f);
dfdx = zeros(n, 1);

dfdx(1) = (f(2) - f(1)) / (x(2) - x(1));
for ii = 2:n-1
    dfdx(ii) = (f(ii+1) - f(ii-1)) / (x(ii+1) - x(ii-1));
end
dfdx(n) = (f(n) - f(n-1)) / (x(n) - x(n-1));
end

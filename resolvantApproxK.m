function [ out ] = resolvantApproxK( K, M, b, tau, eta, nb_coeff )
%resolvantApprox Computes a chebyshev approx of (inv(M)*K-tau*I)^{-1} * b
%   Input :
%   - K : symmetric positive matrix
%   - M : symmetric positive definite matrix
%   - b : rhs
%   - tau : complex parameter
%   - xi : real negative shift
%   - nb_coeff : degree of the chebyshev approximant
%   Output :
%   - out : approximation of (inv(M)*K-tau*I)^{-1} * b.

[n ~] = size(K);
z = 2*eta/(eta-tau)-1;

w = z + sqrt(z^2 - 1);   % |w| is the inverse convergence factor
rac = sqrt(z^2 - 1);
if abs(w) < 1,
    w = z - sqrt(z^2 - 1);   % |w| is the inverse convergence factor
    rac = - sqrt(z^2 - 1);
end
an = 2./(rac * w.^(0:nb_coeff));
an(1) = an(1)/2;

b = (K-eta*M) \ (M*b);

T0 = b;
T1 = -2*eta*((K-eta*M)\(M*b)) - b;
out = an(1)*T0 + an(2)*T1;
for j = 3:length(an),
    [T0, T1] = deal(T1, -4*eta*((K-eta*M)\(M*T1)) -2*T1 - T0);
    out = out + an(j)*T1;
end

out = ((2*eta)/(eta-tau)) * out;
end

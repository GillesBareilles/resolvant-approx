function [ out ] = resolvantApprox(K, M, b, xis, tols, cheb_max_order)
%resolvantApprox Computes approximation to (inv(M)*K-xis_i*I)^{-1} * b to a
%tols(i) error.
%   Input :
%   - K : symmetric positive matrix
%   - M : symmetric positive definite matrix
%   - b : rhs
%   - xis : complex shifts
%   - tols : error target for each system
%   - cheb_max_order : maximum degree of Chebyshev approximant allowed
%   Output :
%   - out : approximation of (xis_i*I-inv(M)*K)^{-1} * b.

if nargin ~= 6
    error('Wrong number of input arguments.\n resolvantApprox(K, M, b, xis, tols, cheb_max_order)')
end

if (size(K, 1)~=size(K, 1)) || (size(M, 1)~=size(M, 1))
    error('Input matrices must be square')
end
if (size(K)~=size(M))
    error('Input matrices must have same size')
end
if (size(K, 1)~=size(b, 1) || (size(b, 2)~=1))
    error('Input b must be a vector compatible with input matrices')
end

if ~((size(xis, 1)==1) || (size(xis, 2)==1))
    error('Input xis must be a row or column vector')
end
if ~((size(tols, 1)==1) || (size(tols, 2)==1))
    error('Input tols must be a row or column vector')
end

if size(xis, 1)==1
    xis = transpose(xis);
end
if size(tols, 1)==1
    tols = transpose(tols);
end

x0 = min(real(xis));
objectivefun = @(x) compute_cheb_maxdegree(x, xis, tols);

[eta,~,~,~] = fminsearch(objectivefun,x0);
%fprintf('Minimum found at %f, value %f\n', eta, fval);

cheb_degrees = zeros(length(xis), 1);
for i=1:length(xis)
    cheb_degrees(i) = ceil(compute_cheb_degree(eta, xis(i), tols(i)));
end

for i=1:length(cheb_degrees)
    if(cheb_degrees(i) > cheb_max_order)
        warning('Pole n°%i: required degree for tol %1.1e is %i, brought to %i.', i, tols(i), cheb_degrees(i), cheb_max_order);
        cheb_degrees(i) = cheb_max_order;
    end
end

out = chebyshevEvaluation( K, M, b, xis, eta, cheb_degrees );
end


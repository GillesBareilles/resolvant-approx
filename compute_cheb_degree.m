function [ m_max ] = compute_cheb_degree( eta, tau, tol )
%compute_cheb_degree Computes the minimum chebyshev degree required to achieve the required tolerance for the given xi and eta
%   Input :
%       - eta : real negative number
%       - xi : complex pole
%       - tol : positive number

z = 2*eta/(eta-tau) - 1;
temp_dc = z^2 - 1;

w = z + sqrt(temp_dc); % 1/|w| convergence factor
if (abs(w) < 1)
    w = z - sqrt(temp_dc);
end
m_max = log(4 / (tol*(abs(w)-1)*sqrt(abs(z*z-1))*abs(eta-tau))) / log(abs(w));
end

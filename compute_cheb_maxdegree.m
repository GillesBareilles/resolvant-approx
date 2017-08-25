function [ m_max ] = compute_cheb_maxdegree( eta, taus, tols )
%compute_cheb_maxdegree Computes the maximum chebyshev degree required to achieve the required tolerance for the given xis and eta
%   Input :
%       - eta : real negative number
%       - xis : sequence of complex poles
%       - tols : sequence of real positive numbers

ms = zeros(length(taus), 1);
for i=1:length(taus)
    z = 2*eta/(eta-taus(i)) - 1;
    temp_dc = z^2 - 1;
    
    w = z + sqrt(temp_dc); % 1/|w| convergence factor
    if (abs(w) < 1)
        w = z - sqrt(temp_dc);
    end
    ms(i) = log(4 / (tols(i)*(abs(w)-1)*sqrt(abs(z*z-1))*abs(eta-taus(i)))) / log(abs(w));
end
m_max = max(ms(:));
end


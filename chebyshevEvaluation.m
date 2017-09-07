function [ out ] = chebyshevEvaluation( K, M, b, taus, eta, cheb_degrees )
%chebyshevEvaluation Computes a chebyshev approx of (tau*I-inv(M)*K)^{-1} * b
%   Input :
%   - K : symmetric positive matrix
%   - M : symmetric positive definite matrix
%   - b : rhs
%   - taus : complex parameters
%   - eta : real negative shift
%   - cheb_degrees : degrees of the chebyshev approximants
%   Output :
%   - out : approximation of (tau*I-inv(M)*K)^{-1} * b.

alphas = zeros(size(taus, 2), max(cheb_degrees)+1);
for i=1:length(taus)
    tau = taus(i);
    
    z = 2*eta/(eta-tau)-1;
    
    w = z + sqrt(z^2 - 1);   % |w| is the inverse convergence factor
    rac = sqrt(z^2 - 1);
    if abs(w) < 1,
        w = z - sqrt(z^2 - 1);
        rac = - sqrt(z^2 - 1);
    end
    alphas(i,:) = 2./(rac * w.^(0:(max(cheb_degrees))));
    alphas(i,1) = alphas(i,1)/2;
    % opposite of initial (inv(M)*K-tau*I)^{-1} * b and factor depending on
    % the complex shift
    alphas(i,:) = - alphas(i,:) * ((2*eta)/(eta-tau)); 
end

out = chebyProj(K, M, b, eta, alphas, cheb_degrees);
end

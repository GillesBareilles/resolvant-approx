clear all
close all

compile 

J = 15; h = 2/J;
K = 0.02*gallery('poisson',J-1)/h^2;      % *negative* 2D Laplacian (positive spectrum)
n = length(K);
M = 3*speye(n, n);

%% Testing chebyProj:
fprintf('Testing mex routine...\n');

%Test 1. several "systems", vector rhs
b = rand(n, 1); b = b/norm(b);
alphas = rand(3,20) + 1i*rand(3,20);
cheb_orders = ones(size(alphas, 1), 1) * size(alphas, 2);
eta = -1;

out_c = chebyProj(K, M, b, eta, alphas, cheb_orders);

out = zeros(n, length(cheb_orders));
b = (K-eta*M) \ (M*b);
T0 = b;
T1 = -2*eta*((K-eta*M)\(M*b)) - b;
for i=1:length(cheb_orders)
    out(:,i) = alphas(i,1)*T0 + alphas(i,2)*T1;
end
for j = 3:length(alphas),
    [T0, T1] = deal(T1, -4*eta*((K-eta*M)\(M*T1)) -2*T1 - T0);
    for i=1:length(cheb_orders)
        out(:,i) = out(:,i) + alphas(i,j)*T1;
    end
end

if norm(out_c - out) < 1e-13
    fprintf('\t test 1. ok\n');
else
    fprintf('\t test 1. failed <----\n');
end

%Test 2. one "system", matrix right-hand side b
b = rand(n, 5); b = b/norm(b);
alphas = rand(1,20) + 1i*rand(1,20);
cheb_orders = ones(size(alphas, 1), 1) * size(alphas, 2);
eta = -1;

out_c = chebyProj(K, M, b, eta, alphas, cheb_orders);

out = zeros(n, 1);
b = (K-eta*M) \ (M*b);
T0 = b;
T1 = -2*eta*((K-eta*M)\(M*b)) - b;
out = alphas(1)*T0 + alphas(2)*T1;
for j = 3:length(alphas),
    [T0, T1] = deal(T1, -4*eta*((K-eta*M)\(M*T1)) -2*T1 - T0);
    out = out + alphas(j)*T1;
end

if norm(out_c - out) < 1e-13
    fprintf('\t test 2. ok\n');
else
    fprintf('\t test 2. failed <----\n');
end

%% Testing chebyshevEvaluation:
fprintf('Testing chebyshevEvaluation...\n');

% Simplest test: one pole, one rhs (right hand side b)
b = rand(n, 1); b = b/norm(b);
xis = [-1+0.5i];
eta = -1;
cheb_degrees = 25;

error = norm(chebyshevEvaluation( K, M, b, xis, eta, cheb_degrees ) - (xis*speye(n,n)-inv(M)*K) \ b);

if error < 1e-15
    fprintf('\t test 1. ok\n');
else
    fprintf('\t test 1. failed <----\n');
end

% Several poles, one rhs
b = rand(n, 1); b = b/norm(b);
xis = [-1+0.5i, -2+0.5i];
eta = -1;
cheb_degrees = [25; 25];

out = chebyshevEvaluation( K, M, b, xis, eta, cheb_degrees );
out_m = zeros(n, length(xis));
for i=1:length(xis)
    out_m(:,i) = (xis(i)*speye(n,n)-inv(M)*K) \ b;
end
error = norm(out - out_m);

if error < 1e-15
    fprintf('\t test 2. ok\n');
else
    fprintf('\t test 2. failed <----\n');
end

% Several poles, several rhs
b = rand(n, 5); b = b/norm(b);
xis = [-1+0.5i, -2+0.5i];
eta = -1;
cheb_degrees = [25; 25];

out = chebyshevEvaluation( K, M, b, xis, eta, cheb_degrees );
out_m = zeros(n, length(xis)*size(b, 2));
for i=1:length(xis)
    out_m(:,(size(b, 2)*(i-1)+1):(size(b, 2)*(i))) = (xis(i)*speye(n,n)-inv(M)*K) \ b;
end
error = norm(out- out_m);

if error < 1e-15
    fprintf('\t test 3. ok\n');
else
    fprintf('\t test 3. failed <----\n');
end

%% Testing resolvantApprox
fprintf('Testing chebyshevEvaluation...\n');

% Testing for one rhs, four poles, uniform precision
b = rand(n, 1); b = b/norm(b);
xis = [-1+0.5i, -2+0.5i, -4-0.5i, -1.5-1i];
tols = ones(length(xis), 1) * 1e-9;
cheb_max_order = 25;

out = resolvantApprox(K, M, b, xis, tols, cheb_max_order);

out_m = zeros(length(b), length(xis));
for i=1:length(xis)
    out_m(:,i) = (xis(i)*speye(n,n)-inv(M)*K) \ b;
end

pass = 1;
for i=1:length(xis)
    pass = pass & (norm(out(:,i) - out_m(:,i)) < tols(i));
end

if pass
    fprintf('\t test 1. ok\n');
else
    fprintf('\t test 1. failed <----\n');
end

% Testing for four poles, multiple precisions
b = rand(n, 1); b = b/norm(b);
xis = [  -1.5251e+02 - 5.7042e+02i  -1.5251e+02 + 5.7042e+02i   -1.0653e+02 - 1.1069e+02i   -1.0653e+02 + 1.1069e+02i   -3.6358e+01 - 1.4335e+01i   -3.6358e+01 + 1.4335e+01i   -5.7416e+00 - 4.0420e+00i   -5.7416e+00 + 4.0420e+00i   -6.6606e+02 - 1.8784e+03i   -6.6606e+02 + 1.8784e+03i   -1.1733e+01 + 0.0000e+00i   -8.6931e+00 + 0.0000e+00i ];
tols = ones(12, 1) * 1e-6;
tols(4:6) = 1e-8; tols(7:9) = 1e-10; tols(10:12) = 1e-12;
cheb_max_order = 55; %max required degree should be 52

out = resolvantApprox(K, M, b, xis, tols, cheb_max_order);

out_m = zeros(length(b), length(xis));
for i=1:length(xis)
    out_m(:,i) = (xis(i)*speye(n,n)-inv(M)*K) \ b;
end

pass = 1;
for i=1:length(xis)
    pass = pass & (norm(out(:,i) - out_m(:,i)) < tols(i));
end

if pass
    fprintf('\t test 2. ok\n');
else
    fprintf('\t test 2. failed <----\n');
end

% Testing for four poles, multiple precisions, several poles
b = rand(n, 10); b = b/norm(b);
xis = [  -1.5251e+02 - 5.7042e+02i  -1.5251e+02 + 5.7042e+02i   -1.0653e+02 - 1.1069e+02i   -1.0653e+02 + 1.1069e+02i   -3.6358e+01 - 1.4335e+01i   -3.6358e+01 + 1.4335e+01i   -5.7416e+00 - 4.0420e+00i   -5.7416e+00 + 4.0420e+00i   -6.6606e+02 - 1.8784e+03i   -6.6606e+02 + 1.8784e+03i   -1.1733e+01 + 0.0000e+00i   -8.6931e+00 + 0.0000e+00i ];
tols = ones(12, 1) * 1e-6;
tols(4:6) = 1e-8; tols(7:9) = 1e-10; tols(10:12) = 1e-12;
cheb_max_order = 52; %max required degree should be 52

out = resolvantApprox(K, M, b, xis, tols, cheb_max_order);

out_m = zeros(length(b), length(xis));
for i=1:length(xis)
    out_m(:,(size(b, 2)*(i-1)+1):(size(b, 2)*(i))) = (xis(i)*speye(n,n)-inv(M)*K) \ b;
end

pass = 1;
for i=1:length(xis)
    error = norm(out(:,(size(b, 2)*(i-1)+1):(size(b, 2)*(i))) - (xis(i)*speye(n,n)-inv(M)*K) \ b);
    pass = pass & (error < tols(i));
end

if pass
    fprintf('\t test 3. ok\n');
else
    fprintf('\t test 3. failed <----\n');
end

%% graph ?
% etas = -10:0.1:-0.2;
% ms = zeros(length(etas), 1);
% for i=1:length(xis)
%     for j=1:length(etas)
%         ms(j) = compute_cheb_degree(etas(j), xis(i), tols(i));
%     end
%     plot(etas, ms, '--');
%     hold on
% end
% for j=1:length(etas)
%     ms(j) = compute_cheb_maxdegree(etas(j), xis, tols);
% end
% plot(etas, ms, 'r');
% 
% 
% title('minimum m_i(\eta) for achieving tol_i')
% xlabel('\eta')

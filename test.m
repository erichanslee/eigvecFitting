% Tests fitEigvec routine for stability


% Initialize Data
n = 100; % Size of Matrix
k = 10; % Size of Subset
A = rand(n);
[V, D] = eig(A);
d = diag(D);
win = 1:k;

% ~~~~~ FIRST WORK ON UNSTABLE METHOD ~~~~~~~ %
method = 'both';

% Fit without noise
V_in = V(1:10, win);
V_in = normalizematrix(V_in);
d_in = d(win);
[fRes, fVec] = fitEigvec(A, method, d_in, V_in, win);
display('Clean Eigenvalues:')
disp(d_in)

% Fit with noise
noise = 0.01;
VN_in = V_in + 0.001*randn(size(V_in)).*V_in;
VN_in = normalizematrix(VN_in);
dN_in = d_in + 0.01*randn(size(d_in)).*d_in;
[fResN, fVecN] = fitEigvec(A, method, dN_in, VN_in, win);
display('Noisy Eigenvalues:')
disp(dN_in)

display('Size of non-noisy residual:')
disp(norm(fRes));
display('Size of noisy residual:')
disp(norm(fResN));
display('Quality of Fit of x2:')
disp(abs(diag(normalizematrix(fVec)'*normalizematrix(fVecN))));

% ~~~~~ TRY MORE STABLE METHOD ~~~~~ %
method = 'forward';

% Fit without noise
[fRes, fVec] = fitEigvec(A, method, d_in, V_in, win);
display('Clean Eigenvalues:')
disp(d_in)

% Fit with noise
[fResN, fVecN] = fitEigvec(A, method, dN_in, VN_in, win);
display('Noisy Eigenvalues:')
disp(dN_in)

display('Size of non-noisy residual:')
disp(norm(fRes));
display('Size of noisy residual:')
disp(norm(fResN));

% Things to note:
% -Fitting the usual way leads to far larger error than the newer (but simpler) ways
% -The rest of the eigenvector fit (x2) looks pretty much the same
function [score, multiplier] = stabilityMeasure(noise, method)

% Initialize Data
n = 200; % Size of Matrix
k = 10; % Size of Subset
A = rand(n);
[V, D] = eig(A);
d = diag(D);
win = 1:k;
V_in = V(1:10, win);
V_in = normalizematrix(V_in);
d_in = d(win);
VN_in = V_in + noise*randn(size(V_in)).*V_in;
VN_in = normalizematrix(VN_in);
dN_in = d_in + noise*randn(size(d_in)).*d_in;
switch method
    
    case 'both'
        [fResN, fVecN] = fitEigvec(A, method, dN_in, VN_in, win);
        score = sum(fResN);
        multiplier = max(abs(eigs(A,1,'SM')), abs(eigs(A,1,'LM')));
        
    case 'forward';
        [fResN, fVecN] = fitEigvec(A, method, dN_in, VN_in, win);
        score = sum(fResN);
        multiplier = 1;
        
end


end
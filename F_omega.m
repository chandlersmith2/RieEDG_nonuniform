function X = F_omega(X_in,samples)
% 
    n = size(X_in,1);
    L = n*(n-1)/2;
    I = samples(:,1);
    J = samples(:,2);
    edgeind = I + (J - 1)*n;
    diagind = (1:n + 1:n*n);
    X_diag= diag(X_in);
    Z1 = X_diag(I)+X_diag(J)-2*X_in(edgeind);
    X = zeros(n,n);
    X(edgeind) = -Z1(1:end);
    X(diagind) =  -sum(X,2);
%     X = X./2;

return
   
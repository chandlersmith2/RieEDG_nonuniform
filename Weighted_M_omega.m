function X = Weighted_M_omega(X,samples,p_list)

    n = size(X,1);
    D = diag(X)*ones(1,n) + ones(n,1)*diag(X)' - 2*X;
    P_omega_D = Weighted_P_omega(D,samples,p_list);

    JDJ = -1/2*(P_omega_D - 1/n*ones(n,1)*sum(P_omega_D,2)' - 1/n*sum(P_omega_D,1)'*ones(1,n) + 1/n^2*sum(sum(P_omega_D))*ones(n));
    P_JDJ = Weighted_P_omega(JDJ,samples,p_list);
    RpsRp_X = P_JDJ - diag(sum(P_JDJ,1));

    I = samples(:,1);
    J = samples(:,2);
    edgeind = I + (J - 1)*n;
    diagind = (1:n + 1:n*n);
    X_diag= diag(X);
    Z1 = X_diag(I)+X_diag(J)-2*X(edgeind);
    FpX = zeros(n,n);
    FpX(edgeind) = Z1(1:end).*(p_list.^(-1) - p_list.^(-2));
    FpX(diagind) =  -sum(FpX,2);

    X = RpsRp_X + FpX;
    
return
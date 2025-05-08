function [X_l,P_approx,norm_diffs] = distgeo_rgrad_Mpo(X_0,samples,r,max_iter,M_omega_X,rel_thresh,X,p_list)
%INPUTS:
% X_0: Initialization. In the testing script, this is either the one step
% hard thresholding, or the resampled initialization. Note that the hard
% thresholding step is done in this algorithm, and thus does not need to be
% performed in any pre-processing.
%samples: This is an m by 2 array of sampled indices. It needs to be
%symmetric.
%d: Rank of the manifold for the descent procedure
%max_iter: Max number of iterations before the algorithm is terminated
%F_omega_X: This is the visible information from the sampling, and is
%important for the gradient
%rel_thresh: This is another stopping condition, this one is for relative
%Frobenius norm difference between iterates.
%X: This is the ground truth, and is not necessary for the algorithm. If
%one wishes to compute the convergence plots, returned as norm_diffs, one
%can pass the ground truth in as X.

%OUTPUTS:
%X_l: Final gram matrix recovered
%P_approx: MDS is performed on the final X_l to give the final point cloud
%embedding in d dimensions
%norm_diffs: This returns the relative difference in frobenius norm of X_l
%and X at each iterate l

norm_diffs = zeros(max_iter,1);
%initialize the iteration using the hard thresholding operator
[X_l,~,D_l,U_l]=hard_thresh(X_0,r,1);

for l = 1:max_iter
    % Construct the Euclidean gradient
    G_l = (M_omega_X - Weighted_M_omega(X_l, samples, p_list));

    % Project to the tangent space
    PU_Gl = U_l * (U_l' * G_l);
    Pt_Gl = PU_Gl + PU_Gl' - (PU_Gl * U_l) * U_l';

    % Line search for step size
    Momega_Pt_Gl = Weighted_M_omega(Pt_Gl, samples, p_list);
    alpha_l = (norm(Pt_Gl, 'fro')^2) / sum(sum(Pt_Gl .* Momega_Pt_Gl));

    % Update
    Z = alpha_l * G_l;
    ZU = Z * U_l;
    Y1 = ZU - (U_l * (U_l' * ZU));
    [Q,R] = qr(Y1, 'econ');
    M_l = [D_l + U_l' * ZU, R'; R, zeros(r)];

    % Hard threshold back to manifold
    [U, diagD] = eig(M_l);
    diagD = real(diag(diagD));
    [diagD, idx] = sort(diagD, 'descend');
    diagD = max(diagD(1:r), 0);
    D_l = diag(diagD);

    % Build new iterate
    U = U(:, idx(1:r));
    U_l = [U_l Q] * U;
    X_l1 = X_l;
    X_l = U_l * D_l * U_l';

    norm_diffs(l) = norm(X_l - X, 'fro') / norm(X, 'fro');
    if norm(X_l1 - X_l, 'fro') / norm(X_l, 'fro') < rel_thresh
        break
    end
end

%Perform MDS to return the point cloud representation
P_approx = classical_edg(X_l,r,0);

return


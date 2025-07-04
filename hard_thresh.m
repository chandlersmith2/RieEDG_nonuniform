function [X,U,D,V] = hard_thresh(Y,r,mode)

    % Use partial decompositions to avoid full matrix factorizations.

if mode == 1
    % For symmetric matrices, use eigs for partial eigen-decomposition
    [V,D] = eigs(Y,r,'largestreal');
    X = V*D*V';
    U = 0;
else
    % For general matrices, use svds for partial SVD
    [U,S,V] = svds(Y,r);
    X = U*S*V';
    D = S;
end
return

% % 1 is eig mode, 2 is svd mode, anything else gives full svd

% if mode == 1
%     %for symmetric matrices, do an eigenvalue decomposition in place of an svd
% %     [V,D]=eigs(Y,r,'largestreal');
% %     X = V*D*V';
% %     U = 0;
%     [V,D] = eig(Y);
%     D = diag(D);
%     [~,IJ] = sort(abs(D),'descend');
%     D = D(IJ);
%     V = V(:,IJ);
%     V = V(:,1:r);
%     D = diag(D(1:r));
%     X = V*D*V';
%     U=0;
% else
%     %for general matrices, use the svd instead to preserve the signs of the
%     %singular vectors!
%     [U,D,V]=svd(Y);
%     D = diag(D);
%     D = D(1:r);
%     D = diag(D);
%     U = U(:,1:r);
%     V = V(:,1:r);
%     X = U*D*V';
% end

% return

function Y = Weighted_P_omega(X,samples,p_list)


    n = size(X,1);
    Y = zeros(n);
    IJ = sub2ind([n n],samples(:,1),samples(:,2));
    Y(IJ) = ((p_list).^(-1)).*X(IJ); 

return
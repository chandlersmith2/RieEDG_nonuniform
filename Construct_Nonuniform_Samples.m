function [samples, Weight, p_list] = Construct_Nonuniform_Samples(X,P,rate,distance_mode,lambda)

n = size(X,1);
L = n*(n-1)/2;
m = rate*L;
r = size(P,2);
Weight = rand(n);
[U,E] = eigs(X,r);
nuc_norm_X = sum(abs(diag(E)));
p_list = zeros(n);

if distance_mode == 1
    for i = 1:n
        for j = (i+1):n
            p_ij = (lambda*norm(P(i,:)-P(j,:))^2/((n+1)*nuc_norm_X) + (1-lambda)/L)*m;
            if Weight(i,j) < p_ij
                Weight(i,j) = 1;
                p_list(i,j) = p_ij;
            else
                Weight(i,j) = 0;
            end
        end
    end
elseif distance_mode == 0
    for i = 1:n
        for j = (i+1):n
            p_ij = (lambda*norm(U(i,:)-U(j,:))^2/((n+1)*r) + (1-lambda)/L)*m;
            if Weight(i,j) < p_ij
                Weight(i,j) = 1;
                p_list(i,j) = p_ij;
            else
                Weight(i,j) = 0;
            end
        end
    end
end

Weight = triu(Weight);
Weight = Weight + Weight';
Weight(logical(eye(size(Weight)))) = 0;
[I,J] = find(Weight==1);
samples = [I,J];

p_list = triu(p_list) + triu(p_list,1)'; % Make p_list symmetric
p_list = p_list(sub2ind(size(p_list), samples(:,1), samples(:,2))); % Extract entries corresponding to samples
p_list(p_list > 1) = 1;
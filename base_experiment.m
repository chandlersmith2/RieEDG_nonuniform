% Select the datasets
data_array = ["./synthetic_datasets/cow.off"];%,"./synthetic_datasets/cow.off","swiss","cities"];

% Choose the sampling rates
d_values = [5,10,15,20];
% Number of times to run each dataset (with a new randomized subset of
% points)
num_trials = 1;

% Max iterations stopping condition for the algorithms
max_iter = 1000;

% relative difference threshold 
rel_thresh = 1e-7;
disp_progress = 1;

lsopts.maxit = 20;
lsopts.xtol = 1e-8;
lsopts.gtol = 1e-8; 
lsopts.ftol = 1e-10; 
lsopts.alpha  = 1e-3; 
lsopts.rho  = 1e-4; 
lsopts.sigma  = 0.1; 
lsopts.eta  = 0.8; 

opts.maxit = ceil(max_iter/lsopts.maxit);
opts.r = 1.0;
opts.printenergy = 0;
opts.printerror = 0;
opts.tol = rel_thresh;


F_IPM_results = zeros(length(data_array),length(d_values),num_trials);
F_RMSE_results = zeros(length(data_array),length(d_values),num_trials);
ac_IPM_results = zeros(length(data_array),length(d_values),num_trials);
ac_RMSE_results = zeros(length(data_array),length(d_values),num_trials);

tiledlayout(1,4)

for i=1:length(data_array)
    data_loc = data_array(i);
    if data_loc == "./synthetic_datasets/cow.off" %datapoints P and mesh trg
        [P, trg] = ReadOFF(data_loc,'1');
    elseif data_loc == "./synthetic_datasets/1k.off"
        [P, trg] = ReadOFF(data_loc,'1'); %datapoints P and mesh trg
    elseif data_loc == "swiss"
        load('./synthetic_datasets/ptswiss.mat');
        P = pt;
    elseif data_loc == "cities"
        load('./synthetic_datasets/UScities.mat');
        P = spt(:,1:2); %can make the last arg : to get altitude component
    else
        disp 'You should write something here to load your dataset!'
        break
    end
    
    %number of datapoints
    n = size(P,1);
    %dimension of datapoints
    rank = size(P,2);
    opts.rank = rank;
    %number of squared distances (useful for one of the sampling modes)
    L = n*(n-1)/2;
    %make sure that the datapoints have zero mean for reconstruction
    P = P - sum(P,1)/n; %centering the data

    %build the true gram and distance matrices
    X = P*P';%gram matrix
    D = ones(n,1)*diag(X)'+diag(X)*ones(1,n)-2*X;%distance matrix
    for q=1:n
        D(q,q) = 0.0;
    end

    %cycle through rates
    for j=1:length(d_values)
        d = d_values(j);
        for k=1:num_trials
            Weight = createRandRegGraph(n, d);
            [I,J] = find(Weight==1);
            F_omega_X = F_omega(X,[I,J]);
            X_0 = R_omega(X,[I,J]);

%             [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Fomega(X_0, ...
%                                             [I,J],rank,max_iter, ...
%                                             F_omega_X,rel_thresh,X); 
%             IPM_err = norm(X_approx-X,'fro')/norm(X,'fro');
%             F_IPM_results(i,j,k) = IPM_err; 
%                     
%             P_edg = P_approx(:,1:rank);
%             P_edg = P_edg - mean(P_edg);
%             [U,~,V] = svd(P_edg'*P);
%             RR = U*V';
%             P_edg_rot = P_edg*RR;
%             RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));
%             F_RMSE_results(i,j,k) = RMSE;
% 
%             if disp_progress == 1
%                 disp("F_omega " + data_loc + " d value " + string(d) + ...
%                     " "+string(k)+"th trial IPM relative difference = " ...
%                     + string(IPM_err))
%                 disp(string(k)+"th trial RMSE = " + string(RMSE))
%             end

            [~, X_approx, ~] = alternating_completion(D,Weight,opts,lsopts);
            IPM_err = norm(X - X_approx, 'fro') / norm(X, 'fro');
            [V2,Lam2] = eigs(X_approx,rank,'lm');
            lam2 = diag(Lam2);
            [lam2,IJ2] = sort(lam2,'descend');
            V2 = V2(:,IJ2);
            Pt_edg = real(V2(:,1:rank)*diag(sqrt(lam2(1:rank))));
            pt2 = Pt_edg - mean(Pt_edg);
            [U,~,V] = svd(pt2'*P);
            RR = U*V';
            pt2 = pt2*RR;
            
            RMSE = sqrt(mean(sum((pt2 - P).^2,2)));

            ac_IPM_results(i,j,k) = IPM_err; 
            ac_RMSE_results(i,j,k) = RMSE; 

            if disp_progress == 1
                disp("AC " + data_loc + " d value " + string(d) + ...
                    " "+string(k)+"th trial IPM relative difference = " ...
                    + string(IPM_err))
                disp(string(k)+"th trial RMSE = " + string(RMSE))
            end

            nexttile
            scatter(pt2(:,1),pt2(:,2),'k.')
            
        end
    end
end




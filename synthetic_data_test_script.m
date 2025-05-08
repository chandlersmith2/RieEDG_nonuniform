
% This determines which algorithm you're using. Set it to
% 1 for F_omega, 2 for R_omega
descent_mode = 2;

%1 for one step hard thresholding, 2 for Riemannian Resampling Initialization
initialization_type = 1;

%set sampling_mode = 0 for Bernoulli sampling of the entries, sampling_mode
%= 1 for a uniform with replacement sampling model
sampling_mode = 1;

%select 1 if you want to print progress at each iteration
disp_progress = 1;

% Select the datasets
data_array = ["./synthetic_data/cow.off","./synthetic_data/cow.off","swiss","cities"];

% Choose the sampling rates
rate_array = [.1,.07,.05,.03,.02,.01];

%if using the Resampling initialization, give an estimate for your
%coherence parameters and set your partition size
coherence = 50;
partition_size = 3;

% Number of times to run each dataset (with a new randomized subset of
% points)
num_trials = 2;

% Max iterations stopping condition for the algorithms
num_iter = 50;

% relative difference threshold 
rel_thresh = 10^-7;

% 3d result array
IPM_results = zeros(length(data_array),length(rate_array),num_trials);
RMSE_results = zeros(length(data_array),length(rate_array),num_trials);

for i=1:length(data_array)
    data_loc = data_array(i);
    if data_loc == "./data/cow.off" %datapoints P and mesh trg
        [P, trg] = ReadOFF(data_loc,'1');
    elseif data_loc == "./data/1k.off"
        [P, trg] = ReadOFF(data_loc,'1'); %datapoints P and mesh trg
    elseif data_loc == "swiss"
        load('./data/ptswiss.mat');
        P = pt;
    elseif data_loc == "cities"
        load('./data/UScities.mat');
        P = spt(:,1:2); %can make the last arg : to get altitude component
    else
        disp 'You should write something here to load your dataset!'
        break
    end
    
    %number of datapoints
    n = size(P,1);
    %dimension of datapoints
    d = size(P,2);
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

    figure
    tiledlayout(length(rate_array),num_trials) 
    title(data_loc)
    %cycle through rates
    for j=1:length(rate_array)
        rate = rate_array(j);
        for k=1:num_trials       

            %build the random samples of the matrix
            %samples is the list passed through to the algorithm, Weight is
            %a binary masking matrix that indicates where entries have been
            %sampled. Use spy(Weight) to see where your entries are
            %plotted!
            if initialization_type ~= 2
                [samples,Weight] = Construct_Samples(n,n,L,rate,sampling_mode,1);
            else
                [samples,~] = Construct_Samples(n,n,L,rate,sampling_mode,0);
            end

            if descent_mode == 1
                if initialization_type == 1
                    %One step hard thresholding
                    sym_samples = samples;
                    m = length(samples);
                    R_omega_X = R_omega(X,sym_samples);
                    X_0 = L/m*R_omega_X;
                    F_omega_X = F_omega(X,sym_samples);

                elseif initialization_type == 2
                    %Riemannian resampling
                    X_0 = riemannian_resampling(X,partition_size,samples, ...
                        coherence,d);
                    I = samples(:,1);
                    J = samples(:,2);
                    A = sub2ind([n n],I,J);
                    Weight = zeros(n);
                    Weight(A) = 1;
                    for l=1:n
                        %diagonal entries are known i.e. D_ii = 0
                        Weight(l,l)= 1;
                        for p=l+1:n
                          %make the weight matrix symmetric
                          Weight(p,l)= Weight(l,p);
                        end
                    end
                    [I,J] = find(Weight==1);
                    sym_samples = [I,J];

                    F_omega_X = F_omega(X,sym_samples);
                else
                    disp 'Only hard thresholding and Riemannian resampling are currently implemented'
                    break
                end
               
                [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Fomega(X_0, ...
                                                sym_samples,d,num_iter, ...
                                                F_omega_X,rel_thresh,X); 
                IPM_err = norm(X_approx-X,'fro')/norm(X,'fro');
                IPM_results(i,j,k) = IPM_err; 
                
                P_edg = P_approx(:,1:d);
                P_edg = P_edg - mean(P_edg);
                [U,~,V] = svd(P_edg'*P);
                RR = U*V';
                P_edg_rot = P_edg*RR;
                RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));
                RMSE_results(i,j,k) = RMSE;

                if disp_progress == 1
                    disp("F_omega " + data_loc + "rate " + string(rate) + ...
                        " "+string(k)+"th trial IPM relative difference = " ...
                        + string(IPM_err))
                    disp(string(k)+"th trial RMSE = " + string(RMSE))
                end
                
                norm_diffs = norm_diffs(norm_diffs>0);
                nexttile
                semilogy(1:length(norm_diffs),norm_diffs)
                xlabel("num iter")
                ylabel("rate = "+string(rate))

            elseif descent_mode == 2
                if initialization_type == 1
                    sym_samples = samples;
                    m = length(samples);
                    R_omega_X = R_omega(X,samples);
                    X_0 = L/m*R_omega_X;
                elseif initialization_type == 2
                    X_0 = riemannian_resampling(X,partition_size,samples,coherence,d);
                    I = samples(:,1);
                    J = samples(:,2);
                    A = sub2ind([n n],I,J);
                    Weight = zeros(n);
                    Weight(A) = 1;
                    for l=1:n
                        %diagonal entries are known i.e. D_ii = 0
                        Weight(l,l)= 1;
                        for p=l+1:n
                          %make the weight matrix symmetric
                          Weight(p,l)= Weight(l,p);
                        end
                    end
                    [I,J] = find(Weight==1);
                    sym_samples = [I,J];
                    R_omega_X = R_omega(X,sym_samples);
                else
                    disp 'Choose either 0 or 1 for initialization'
                    break
                end
                
                [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Romega(X_0, ...
                                                            sym_samples,d,num_iter, ...
                                                            R_omega_X,rel_thresh,X);    

                IPM_err = norm(X_approx-X,'fro')/norm(X,'fro');
                IPM_results(i,j,k) = IPM_err;

                P_edg = P_approx(:,1:d);
                P_edg = P_edg - mean(P_edg);
                [U,~,V] = svd(P_edg'*P);
                RR = U*V';
                P_edg_rot = P_edg*RR;
                RMSE = sqrt(mean(sum((P_edg_rot - P).^2,2)));

                if disp_progress == 1
                    disp("R_omega " + data_loc + "rate " + string(rate) + " " ...
                        +string(k)+"th trial IPM relative difference = " ...
                        + string(IPM_err))
                    disp(string(k)+"th trial RMSE = " + string(RMSE))
                end

                norm_diffs = norm_diffs(norm_diffs>0);
                nexttile
                semilogy(1:length(norm_diffs),norm_diffs)
                xlabel("num iter")
                ylabel("rate = "+string(rate))
                title(data_loc+" trial " + string(k))
            else
                disp('Only two algorithms to choose from!')
            end
        end
    end
end

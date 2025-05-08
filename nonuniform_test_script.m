
%= 1 for a uniform with replacement sampling model
distance_mode = 0;

%select 1 if you want to print progress at each iteration
disp_progress = 1;

% Select the datasets
data_array = ["ellipse"];%,"./data/1k.off"];%,"./data/cow.off","swiss","cities"];

% Choose the sampling rates
rate_array = [.1];%,];

% Number of times to run each dataset (with a new randomized subset of
% points)
num_trials = 1;

% Max iterations stopping condition for the algorithms

% relative difference threshold 
rel_thresh = 10^-10;

max_iter = 500;

lambda = 0.8;

% 3d result array
IPM_results = zeros(length(data_array),length(rate_array),num_trials);

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
    elseif data_loc == "ellipse"
        n = 1000;
        theta = 2 * pi * rand(n, 1);
        phi = pi * rand(n, 1);
        x = cos(theta) .* sin(phi);
        y = sin(theta) .* sin(phi);
        z = 0.5 * cos(phi); % Highly degenerate along the z-axis
        P_ellipse = [x, y, z];
        % Add two outlier points along the z-axis
        outlier1 = [0, 0, 10];
        outlier2 = [0, 0, -10];
        P = [P_ellipse; outlier1; outlier2];
        % P = P_ellipse;
    else
        disp 'You should write something here to load your dataset!'
        break
    end
    
    %number of datapoints
    n = size(P,1);
    %dimension of datapoints
    r = size(P,2);
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
    sgtitle(['Dataset: ', char(data_loc), ', Lambda: ', num2str(lambda)])
    tiledlayout(length(rate_array),num_trials) 
    figure
    sgtitle(['Additional Visualization for Dataset: ', char(data_loc)])
    tiledlayout(length(rate_array), num_trials)
    %cycle through rates
    for j=1:length(rate_array)
        rate = rate_array(j);
        for k=1:num_trials       
            [samples,Weight,p_list] = Construct_Nonuniform_Samples(X,P,rate,distance_mode,lambda);
            Mpo_X = Weighted_M_omega(X,samples,p_list);
            X_0 = 1/rate*R_omega(X,samples);

            [X_approx,P_approx,norm_diffs] = distgeo_rgrad_Mpo(X_0,samples,r,max_iter,Mpo_X,rel_thresh,X,p_list);
            ipm_err = norm(X-X_approx,'fro')/norm(X,'fro');
            IPM_results(i,j,k) = ipm_err;

            P_approx = P_approx - mean(P_approx);
            [U,~,V] = svd(P_approx'*P);
            RR = U*V';
            P_approx = P_approx*RR;

            %plot the results
            figure(1) % Ensure plotting on the first tiledlayout figure
            nexttile
            semilogy(norm_diffs,'LineWidth',2)
            title(['Rate = ',num2str(rate),', Trial = ',num2str(k)])
            xlabel('Iteration')
            ylabel('Relative Difference')
            grid on
            if disp_progress
                disp('Dataset: '+ string(data_loc) +  ', Rate: ' + num2str(rate) + ', Trial: ' + num2str(k) + ', Err: ' + num2str(ipm_err))
            end
            figure(2) % Ensure plotting on the second tiledlayout figure
            nexttile
            scatter3(P_approx(:,1), P_approx(:,2), P_approx(:,3), 10, 'blue', 'filled');
            hold on;
            scatter3(P(:,1), P(:,2), P(:,3), 10, 'red', 'filled');
            hold off;
            title(['Rate = ', num2str(rate), ', Trial = ', num2str(k)])
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            legend('P_{approx}', 'P', 'Location', 'best')
            grid on
        end
    end
end


% n = 1000;
% theta = 2 * pi * rand(n, 1);
% phi = pi * rand(n, 1);
% 
% x = cos(theta) .* sin(phi);
% y = sin(theta) .* sin(phi);
% z = 0.5 * cos(phi); % Highly degenerate along the z-axis
% P_ellipse = [x, y, z];
% 
% % Add two outlier points along the z-axis
% outlier1 = [0, 0, 10];
% outlier2 = [0, 0, -10];
% P = [P_ellipse; outlier1; outlier2];
% 
% % Plot side by side
% figure;
% 
% subplot(1, 2, 1);
% scatter3(P_ellipse(:,1), P_ellipse(:,2), P_ellipse(:,3), 10, 'filled');
% title('Incoherent Point Cloud', 'FontSize', 24);  % Updated title
% xlabel('x'); ylabel('y'); zlabel('z');
% axis square;
% grid on;
% 
% subplot(1, 2, 2);
% scatter3(P(:,1), P(:,2), P(:,3), 10, 'filled');
% hold on;
% scatter3([0 0], [0 0], [10 -10], 50, 'r', 'filled'); % highlight outliers
% title('Coherent Point Cloud', 'FontSize', 24);  % Updated title
% xlabel('x'); ylabel('y'); zlabel('z');
% axis square;
% grid on;
% exportgraphics(gcf, 'outlier_figure.png', 'Resolution', 300);


% Assume P, P_1, and P_2 are already defined as n-by-3 matrices

figure;

% Uniform Sampling Comparison
subplot(1,2,1)
scatter3(P(:,1), P(:,2), P(:,3), 10, 'filled', 'MarkerFaceColor', [0.6 0.6 0.6])
hold on
scatter3(P_1(:,1), P_1(:,2), P_1(:,3), 10, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8])
title('Uniform Sampling', 'FontSize', 20)
axis square
grid on
legend('Original', 'Reconstruction', 'Location','northeast','FontSize', 12)
view(3)

% Coherence Based Sampling Comparison
subplot(1,2,2)
scatter3(P(:,1), P(:,2), P(:,3), 10, 'filled', 'MarkerFaceColor', [0.6 0.6 0.6])
hold on
scatter3(P_approx(:,1), P_approx(:,2), P_approx(:,3), 10, 'filled', 'MarkerFaceColor', [0.8 0.2 0.2])
title('Coherence Based Sampling', 'FontSize', 20)
axis square
grid on
legend('Original', 'Reconstruction', 'Location', 'northeast','FontSize', 12)
view(3)
% 
% % Optional: Save high-resolution figure
% % print('comparison_plot','-dpng','-r300')
% exportgraphics(gcf, 'coherentsampling_vs_uniform.png', 'Resolution', 300)

%            scatter(P_approx(:,1), P_approx(:,2), 10, 'blue', 'filled');
%             hold on;
%             scatter(P(:,1), P(:,2), 10, 'red', 'filled');
%             hold off;
%             title(['Rate = ', num2str(rate), ', Trial = ', num2str(k)])
%             xlabel('X')
%             ylabel('Y')
% 
% figure;
% 
% % Uniform Sampling Comparison
% subplot(1,2,1)
% scatter(P(:,1), P(:,2), 10, 'filled', 'MarkerFaceColor', [0.6 0.6 0.6])
% hold on
% scatter(P_1(:,1), P_1(:,2), 10, 'filled', 'MarkerFaceColor', [0.2 0.4 0.8])
% % title('Uniform Sampling', 'FontSize', 14)
% axis square
% grid on
% legend('Original', 'Reconstruction', 'Location', 'northeast', 'FontSize', 12)
% 
% % Coherence Based Sampling Comparison
% subplot(1,2,2)
% scatter(P(:,1), P(:,2), 10, 'filled', 'MarkerFaceColor', [0.6 0.6 0.6])
% hold on
% scatter(P_2(:,1), P_2(:,2), 10, 'filled', 'MarkerFaceColor', [0.8 0.2 0.2])
% % title('Coherence Based Sampling', 'FontSize', 14)
% axis square
% grid on
% legend('Original', 'Reconstruction', 'Location', 'northeast', 'FontSize', 12)

exportgraphics(gcf, 'coherentsampling_vs_uniform_USCities.png', 'Resolution', 300)
%% This demo file runs EBDR algorithm for p-nearest neighbor selection
% For details, see: 
%
% [1] A. Tastan, M. Muma and A. M. Zoubir, "Eigenvalue-Based Block Diagonal
% Representation and Application to p-Nearest Neighbor Graphs," in Proc.
% 30th European Signal Process. Conf. (accepted), 2022.
%     
%
% Copyright (C) 2022 Aylin Tastan. All rights reserved.
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <https://www.gnu.org/licenses/>.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
clear all
close all

%% Call real-world data set and normalize
[X, num_features, num_samples_total, num_clusters] = call_dataset();
X = normalizing_feature_vectors(X,num_samples_total); %the feature vectors are normalized i.e. ||x||_2

%% Define parameters
TOL = 1e-3;  %Tolerence level
% mix_order = randperm(num_samples_total);
% X = X(:,mix_order);
plotting = 0;
k = num_clusters;
number_cand_max = 8; %number of candidate linear pieces
min_num_samples = num_samples_total/(2*k); %optional for computational cost
MC = 100; %number of Monte Carlo runs


tic
for t = 1:MC
    %% Generate the initial affinity matrix
    [W_ini,dump,idx] = compute_affinity(X); %zero diagonal affinity matrix

    %% p-nearest neighbor selection
    p_set = [5:5:size(X,2)-1]; %Candidate set of p
    
    for i = 1:length(p_set)
        error = inf;

        %Generate the vector of eigenvalues Lambda associated with p_i
        W_p_i = find_affinity_p_nearest_neighbor(dump,idx,p_set(i),num_samples_total); %zero diagonal affinity matrix
        D_p_i = diag(sum(W_p_i)); %Diagonal weight matrix of size n x n
        L_p_i = D_p_i-W_p_i;%Laplacian matrix associated with p_i
        [~,diag_Lambda_p_i] = eig(L_p_i,D_p_i); %eigenvalue decomposition
        [Lambda_p_i,ind_Lambda_p_i] = sort(diag(diag_Lambda_p_i));

        %Compute the vector v associated with p_i
        [v_p_i,L_block_p_i] = compute_vector_v(L_p_i);

        %Compute the matrix of candidate block sizes
        n_cand = compute_candidate_vectors_n(v_p_i,number_cand_max,k,min_num_samples);

        for cand = 1:size(n_cand,1)

            %Estimate the block size vector n associated with p_i
            n_p_i_cand = n_cand(cand,:);

            %Estimate the block coefficient vector w associated with p_i
            [w_p_i_cand,v_hat_p_i_cand] = compute_vector_w(v_p_i,n_p_i_cand,length(n_p_i_cand),plotting);

            %Checking condition
            if(any(w_p_i_cand >= 1)|any(w_p_i_cand <= 0))
                error_cand = inf;
            else
                %Calculate bias
                error_cand = norm(v_p_i-v_hat_p_i_cand)^2;
                clear w_epsilon_i_true
            end

            %Update the error
            if(error_cand < error)
                error = error_cand;
                v_hat_p_i = v_hat_p_i_cand;
                n_p_i = n_p_i_cand;
                w_p_i = w_p_i_cand;
            end
        end

        if (isinf(error))  %Check goodness of the fit
            error_candp(i) = inf;

        else
            %Estimate the ideal vector of eigenvalues 
            Lambda_target_p_i = design_Lambda_target(n_p_i,k);

            %Calculate error
            error_candp(i) = norm(Lambda_p_i-Lambda_target_p_i)^2;
        end


        %Breaking rule
        if (error_candp < TOL )
            W_hat = W_p_i; D_hat = D_p_i; L_hat = L_p_i;
            break;
        end
    end

    %% Estimate p and associated parameters
    %Estimate p-nearest neighbor
    [error_hat,ind_error_hat] = min(error_candp);
    p_hat = p_set(ind_error_hat);

    %Estimate block diagonal affinity matrix
    W_hat = find_affinity_p_nearest_neighbor(dump,idx,p_hat,num_samples_total); %zero diagonal affinity matrix
end
toc
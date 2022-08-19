% This function computes candidate size vectors based on a piece-wise
% linear fit of vector v. For details, see :
%
% [1] A. Tastan, M. Muma and A. M. Zoubir, "Eigenvalue-Based Block Diagonal
% Representation and Application to p-Nearest Neighbor Graphs," in Proc.
% 30th European Signal Process. Conf. (accepted), 2022.
%
% Copyright (C) 2022  Aylin Tastan. All rights reserved.
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
% for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program.  If not, see <https://www.gnu.org/licenses/>.
%
% Inputs          :
% v               : a numeric vector v of size n x n
% number_cand_max : a numeric denoting max number of changepoints
% k               : the number of linear pieces
% min_num_samples : defined minimum number of nodes in blocks
%
% Outputs         :
% n_cand          : a matrix of candidate size vectors
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function N_cand = compute_candidate_vectors_n(v,number_cand_max,k,min_num_samples)

if(nargin < 4)
    min_num_samples = 0;
end

%Find the candidate changepoints in piece-wise linear function
[freq_hat,~] = find(ischange(v,'linear','MaxNumChanges',number_cand_max));
cand_limits = nchoosek(freq_hat,k-1);
ind_cand_changepoints = [ones(size(cand_limits,1),1),cand_limits,(length(v)+1)*ones(size(cand_limits,1),1)];

%Compute matrix of candidate block sizes N_cand
N_cand = diff(ind_cand_changepoints,1,2); % matrix of candidate size vectors

if(min_num_samples)
    [ind_row,~] = find(N_cand<min_num_samples); ind_row=unique(ind_row); % remove the blocks which has smaller than minimum number of samples
    N_cand(ind_row,:) = [];
end

end
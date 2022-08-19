% This function computes the normalized feature vectors ||x||_2. For
% details, see:
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
% Inputs            : 
% X                 :(numeric) data matrix of size m x n 
%                    (m:number of features , n:number of observations)
% num_samples_total : number of observations/samples
% 
% Output            : 
% X                 : (numeric) data matrix of size m x n whose columns are
%                      normalized (m:number of features , n:number of 
%                      observations)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function X = normalizing_feature_vectors(X,num_samples_total)

for j = 1:num_samples_total
    X(:,j) = X(:,j) ./ twonorm(X(:,j));
end

end
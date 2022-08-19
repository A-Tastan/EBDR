% This function computes an affinity matrix. For details, see:
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
% X               : (numeric) data matrix of size m x n
%                   (m: num_features, n: num_samples_total)
% Outputs         :
% W               : (numeric) affinity matrix of size n x n
% dump            : (numeric) a version of W whose rows are sorted in
%                   descending order
% idx             : (numeric) matrix of size n x n whose rows include
%                   indexes associated with the sorted matrix dump
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [W,dump,idx] = compute_affinity(X)

[nFea,nSmp] = size(X);
nBlock = 4000;

for i = 1:ceil(nSmp/nBlock)
    if i == ceil(nSmp/nBlock)
        smpIdx = (i-1)*nBlock+1:nSmp;
        dist = X(:,smpIdx)'*X;
        dist = full(dist);
        [dump idx] = sort(-dist,2); % sort each row

    else
        smpIdx = (i-1)*nBlock+1:i*nBlock;
        dist = X(:,smpIdx)'*X;
        dist = full(dist);
        [dump idx] = sort(-dist,2); % sort each row

    end
end

W = max(dist,dist');
W = W - diag(diag(W)); %zero diagonal affinity matrix
clear dist

end

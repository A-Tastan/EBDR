% This function computes an affinity matrix associated with a p-nearest
% neighbor graph. For details, see:
%
% [1] A. Tastan, M. Muma and A. M. Zoubir, "Eigenvalue-Based Block Diagonal
% Representation and Application to p-Nearest Neighbor Graphs," in Proc.
% 30th European Signal Process. Conf. (accepted), 2022.
%
% [2] D. Cai, X. He and J. Han, "Document clustering using locality 
% preserving indexing" IEEE Trans. Knowl. Data Eng., vol.17, pp. 1624-1637, 
% 2005.
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
% dump              : (numeric) an affinity matrix of size n x n whose rows 
%                     are sorted in descending order
% idx               : (numeric) matrix of size n x n whose rows include
%                     indexes associated with the sorted matrix dump
% p                 : nearest neighbor value
% num_samples_total : total number of samples n
%
% Outputs           :
% W_p               : (numeric) affinity matrix associated with p-nearest
%                     neighbor graph
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function W_p = find_affinity_p_nearest_neighbor(dump,idx,p,num_samples_total)

    nBlock = 4000;
    G = zeros(num_samples_total*(p+1),3);

    for i = 1:ceil(num_samples_total/nBlock)
        if i == ceil(num_samples_total/nBlock)
            smpIdx = (i-1)*nBlock+1:num_samples_total;
            idx = idx(:,1:p+1);
            dump = -dump(:,1:p+1);
            
            G((i-1)*nBlock*(p+1)+1:num_samples_total*(p+1),1) = repmat(smpIdx',[p+1,1]);
            G((i-1)*nBlock*(p+1)+1:num_samples_total*(p+1),2) = idx(:);
            G((i-1)*nBlock*(p+1)+1:num_samples_total*(p+1),3) = dump(:);
        else
            smpIdx = (i-1)*nBlock+1:i*nBlock;
            idx = idx(:,1:p+1);
            dump = -dump(:,1:p+1);
            
            G((i-1)*nBlock*(p+1)+1:i*nBlock*(p+1),1) = repmat(smpIdx',[p+1,1]);
            G((i-1)*nBlock*(p+1)+1:i*nBlock*(p+1),2) = idx(:);
            G((i-1)*nBlock*(p+1)+1:i*nBlock*(p+1),3) = dump(:);
        end
    end

    W_p = sparse(G(:,1),G(:,2),G(:,3),num_samples_total,num_samples_total);
    clear G dist
    W_p = max(W_p,W_p');
    W_p = full(W_p);  W_p = W_p-diag(diag(W_p)); % full zero diagonal affinity matrix
    
end


% This function orders Laplacian matrix using Reverse Cuthill-McKee (RCM)
% algorithm and computes vector v. For details, see :
%
% [1] A. Tastan, M. Muma and A. M. Zoubir, "Eigenvalue-Based Block Diagonal
% Representation and Application to p-Nearest Neighbor Graphs," in Proc.
% 30th European Signal Process. Conf. (accepted), 2022.
%
% [2] E. Cuthill, and J. McKee, "Reducing the bandwidth of sparse symmetric
% matrices," in Proc. 24th Nat. Conf., 1969.
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
% Inputs   :
% L        : Laplacian matrix of size n x n
%
% Outputs  :
% v        : vector v of size n x 1
% L_block  : block diagonal ordered Laplacian matrix of size n x n
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,L_block] = compute_vector_v(L)

% Bandwidth of the initial matrix
[i,j] = find(L);
bw_L = max(i-j) + 1;  % bandwidth of L

% Apply RCM to design a block diagonal matrix
v_per_hat = symrcm(L);

% Permute the Laplacian to separate the blocks
Lp = L(v_per_hat,v_per_hat); % permuted Laplacian

% Sum the permuted laplacian to compute the vector v
[i,j] = find(Lp);
bw_Lp = max(i-j) + 1;  % bandwidth of the permuted Laplacian matrix Lp

% Use the block diagonal matrix with the minimum bandwidth
if(bw_L>bw_Lp)
    v = sum(triu(Lp),2);
    L_block = Lp;
else
    v = sum(triu(L),2);
    L_block = L;
end

end
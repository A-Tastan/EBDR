% This function computes target Lambda vector which is a vector containing
% target eigenvalues in ascending order. For details, see :
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
%
% Inputs        :
% n             : (numeric) block size vector of size k x 1
% k             : number of blocks
%
% Outputs       :
% Lambda_target : the target vector of eigenvalues whose size n x 1
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function Lambda_target = design_Lambda_target(n,k)

%% Designing ideal vector of eigenvalues
% Zeros
Lambda_zeros = zeros(k,1);

% Coefficients
u_l = 0;
Lambda_ndiv = [];

for l = 1:k
    up_bound = n(l)-1;  %n_l-1
    nw(1:up_bound) = n(l)/(n(l)-1); %Generalized eigen-decomposition
    Lambda_ndiv = [Lambda_ndiv,nw];
    nw = [];
end

Lambda_target = [Lambda_zeros.',sort(Lambda_ndiv)];

end
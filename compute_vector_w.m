% This function estimates vector v based on a piece-wise linear fit. The
% changepoints of v define the blocks sizes and the coefficients around
% which the blocks are concentrated. For details, see :
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
% Inputs   :
% v        : (numeric) vector v of size n x 1
% n        : (numeric) block size vector of size k x 1
% k        : number of blocks
% plotting : a binary variable that indicates whether a figure will be
%            generated or not (default is 0)
%
% Outputs  :
% w_hat    : (numeric) estimated similarity coefficients vector of size
%            k x 1
% v_hat    : (numeric) vector v estimate of size n x 1
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [w_hat,v_hat] = compute_vector_w(v,n,k,plotting)

if(nargin < 4)
    plotting = 0;
end

u_l = 0;
for l = 1:k

    % Generate the lth piece
    ell_l = u_l + 1;     u_l = sum(n(1:l));   % lower and upper bounds for the block l
    v_l = v(ell_l:u_l);

    % Form the sample matrix
    % j=ell_l:u_l; %index vector
    j = 1:n(l);
    v_sample_l = [j;v_l.'].'; %samples vector v

    % Compute covariance matrix
    Sigma_l = cov(v_sample_l);    %2x2 matrix

    % Compute mean vector
    Mu_l = mean(v_sample_l).';  %2x1 vector

    % Estimate the normal vector and the bias using eigen-decomposition
    [Eigenvec_Phi,Lambdas_Phi] = eig(Sigma_l);
    [~, ind_Lambdas] = sort(diag(Lambdas_Phi));
    vartheta_l = Eigenvec_Phi(:,ind_Lambdas(1));
    b_l = -vartheta_l.'*Mu_l;

    % Linear regression
    v_l_hat = -(vartheta_l(1)*j+b_l)/vartheta_l(2);

    % Coefficient correspond to the lth block
    w_l_hat = (v_l_hat(end))/(u_l-ell_l);

    % Coefficient vector
    w_hat(l) = w_l_hat;

    % The linear fit
    v_hat(ell_l:u_l) = v_l_hat;

end

if(plotting == 1)
    plot(v);
    hold on
    plot(v_hat);
    hold off
end

end

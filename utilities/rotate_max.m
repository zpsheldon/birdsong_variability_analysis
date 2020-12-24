function V = rotate_max(W)
% V = rotate_max(W)
% 
% MATLAB S7: Takes global factor matrix W and rotates it as described in Glaze &
% Troyer (2012): Maximizes the sum of the values in the first column and
% orthonalizes the remaining columns, sorting by eigenvalue.
%
% AUTHOR: Chris Glaze cmglaze@gmail.com

[K,Q] = size(W);
One = ones(K,1);

% Create first column
V1 = W*W'*One / norm(W'*One);

% Find eigenvectors and eigenvales of residual factor matrix
[Vr,Lambda] = eig(W*W' - V1*V1');

% Sort by eigenvalue, take Q-1 largest, and scale associated eigenvectors
[lambda,sortinds] = sort(diag(Lambda),'descend');
matinds = sortinds(1:Q-1);
Vr = Vr(:,matinds) * sqrt(Lambda(matinds,matinds));

% Concatenate first column with result
V = [V1 Vr];

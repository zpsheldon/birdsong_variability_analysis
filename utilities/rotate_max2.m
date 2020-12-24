function V = rotate_max2(W)
% V = rotate_max(W)
% 
% Takes global factor matrix W and rotates it as described in Glaze &
% Troyer (2012): Maximizes the sum of the values in the first column and
% orthonalizes the remaining columns, sorting by eigenvalue.
%
% AUTHOR: Chris Glaze cmglaze@gmail.com

[K,Q] = size(W);
One = ones(K,1);
Alt = repmat([1;-1],ceil(K/2),1);
Alt = Alt(1:K);

% Create first column
V1 = W*W'*One / norm(W'*One);

% Find eigenvectors and eigenvales of residual factor matrix
[Vr,Lambda] = eig(W*W' - V1*V1');

% Sort by eigenvalue, take Q-1 largest, and scale associated eigenvectors
[lambda,sortinds] = sort(diag(Lambda),'descend');
matinds = sortinds(1:Q-1);
Vr = Vr(:,matinds) * sqrt(Lambda(matinds,matinds));

V2 = Vr*Vr'*Alt / norm(Vr'*Alt);


[Vr,Lambda] = eig(W*W' - V1*V1' - V2*V2');

% Sort by eigenvalue, take Q-1 largest, and scale associated eigenvectors
[lambda,sortinds] = sort(diag(Lambda),'descend');
matinds = sortinds(1:Q-2);
Vr = Vr(:,matinds) * sqrt(Lambda(matinds,matinds));



% Concatenate first column with result
V = [V1 V2 Vr];

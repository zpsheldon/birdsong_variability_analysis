function SRMR = calcSRMR(X,modelstrct)

R = corrcoef(X);
K = length(R);

D = -diff(eye(K))';
Shat = modelstrct.W*modelstrct.W' + D*diag(modelstrct.sigma)*D' + diag(modelstrct.psi);
stdX = std(X);
Rhat = Shat ./ (stdX'*stdX);

% Indexing to unique entries (correlation and covariance matrices are
% symmetric, etc.).
u_inds = find(tril(ones(K)));

SRMR = sqrt(mean((R(u_inds)-Rhat(u_inds)).^2));
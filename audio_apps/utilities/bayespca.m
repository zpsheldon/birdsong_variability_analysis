function [W,sigma,z,alpha,iter,logp] = bayespca(X,m_glob,W0,sigma0)

maxiter = 1000;

if nargin < 2
    m_glob = 1;
end

conv_crit = 0.01;
[n,d] = size(X);

X = X - repmat(mean(X),n,1);

C = X'*X/n;
vr = diag(C);

if nargin < 3
    W0 = repmat(sqrt(vr),1,m_glob).*rand(d,m_glob);
    sigma0 = mean(vr)*rand;
end

W = W0;
sigma = sigma0;

alpha = m_glob ./ diag(W'*W);

if nargout > 6
    logp = zeros(maxiter,1);
end

iter = 1;

sigmamat = zeros(maxiter,d);
wmat = zeros(maxiter,d,m_glob);

sigmamat(iter,:) = sigma;
wmat(iter,:,:) = W;


Chat = W*W' + sigma*eye(d);
logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));


deltaflg = 1;

W_prev = W;
sigma_prev = sigma;

while iter < maxiter && deltaflg

    %     iter

    invsigma = inv(sigma);

    M = W'*W + sigma*eye(m_glob);
    invM = inv(M);

    z = X*W*invM;


    z2sm = z'*z + n*invM*sigma;
    W = (X'*z * inv(z2sm+sigma*diag(alpha)));
    alpha = m_glob ./ diag(W'*W);
    
    sigma = mean(diag(C)) + trace(z2sm*W'*W)/(n*d) - 2*trace(W'*(X'*z))/(n*d);

    iter = iter + 1;

    if nargout > 4
        Chat = W*W' + sigma*eye(d);
        logp(iter) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
    end

    Wdiff = 100 * max(abs(1 - W(:) ./ W_prev(:)));
    
    sigmadiff = 100 * max(abs(1 - sigma ./ sigma_prev(:)));
    
    W_prev = W;
    sigma_prev = sigma;
    
    deltaflg = ~(sigmadiff < conv_crit && Wdiff < conv_crit);


end

if nargout > 6
    logp = logp(1:iter);
end

lat = sum(W.^2);
[dmy,srtinds] = sort(lat,'descend');
W = W(:,srtinds);
z = z(:,srtinds);
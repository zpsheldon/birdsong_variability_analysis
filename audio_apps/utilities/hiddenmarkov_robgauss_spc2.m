function [A,mu,sigma,Pi,v,L,DIC,W,iter,s,Lvc,u,mumu,logu,logdetlambda,mumat,vmat] = hiddenmarkov_robgauss(X,k,alphak,zeroinds,A0,Pi0,mu0,sigma0,v0)

Lconvflg = 0;

[n,d] = size(X);

maxiter = 1000;
C = cov(X);

if nargin < 3
    alphak = 0.001;
end

if nargin < 4
    zeroinds = [];
end

if nargin < 9
    if ~isempty(zeroinds)
        A0 = rand(k);
        A0(1,1) = 0;
        ptmp = sum(A0) / sum(A0(:));
        
        mu0 = repmat(mean(X),k,1) .* ((1/3) + (4/3)*rand(k,d));
        mu0 = repmat(mean(X),k,1) .* mu0 ./ repmat((ptmp*mu0),k,1);
        mu0(1,:) = zeros(1,d);
        
        dmu = mu0 - repmat(mean(X),k,1);
        sigmamu = sum(repmat(ptmp',1,d).*dmu.*dmu);
        
        minvar = diag(diag(diag(C)/(k^2)));
        
        sigma0 = repmat(diag(max(minvar,diag(C)-sigmamu')),[1,1,k]);
        sigma0(:,:,1) = .01*eye(d);
        
        Pi0 = [1 zeros(1,k-1)];
        
    else
        A0 = rand(k);
        A0 = A0 ./ repmat(sum(A0,2),1,k);
        
        ptmp = sum(A0) / sum(A0(:));
        
        mu0 = repmat(mean(X),k,1) .* ((1/3) + (4/3)*rand(k,d));
        mu0 = repmat(mean(X),k,1) .* mu0 ./ repmat((ptmp*mu0),k,1);
        
        dmu = mu0 - repmat(mean(X),k,1);
        sigmamu = sum(repmat(ptmp',1,d).*dmu.*dmu);
        
        minvar = diag(diag(diag(C)/(k^2)));
        
        sigma0 = repmat(diag(max(minvar,diag(C)-sigmamu')),[1,1,k]);
        
        Pi0 = zeros(1,k);
        piconst = (2*pi)^(d/2);
        for kind = 1:k
            Xsub = X(1,:) - mu0(kind,:);
            Pi0(kind) = exp(-.5*diag(Xsub* inv(sigma0(:,:,kind))*Xsub')) / (sqrt(det((sigma0(:,:,kind)))) * piconst);
        end
        
        Pi0 = Pi0 / sum(Pi0);
    end
    
    v0 = (rand(1,k)+1)*20;

    
end

lambda0 = sigma0;
mumu0 = sigma0;

for kind = 1:k
    lambda0(:,:,kind) = inv(sigma0(:,:,kind));
    mumu0(:,:,kind) = mu0(kind,:)'*mu0(kind,:);
end

mu = mu0;
lambda = lambda0;
Pi = Pi0;
mumu = mumu0;
A = A0;

eta0 = 1;
p0 = 0.001;


Winv = zeros(size(lambda));
W = Winv;

alpha0 = alphak * k;

W0 = eye(d);
W0inv = inv(W0);

logdetlambda = zeros(1,k);
v = logdetlambda;

a = zeros(n,k);
b = a;
mandist = a;

% need to initialize logdetlambda, logu, mandist, u, mumu

X1 = repmat(X,[1,1,d]);
X2 = shiftdim(repmat(X',[1,1,d]),1);
Xprod = X1 .* X2;
Xprod = reshape(Xprod,n,d^2);

X2 = reshape(X2,n,d^2);

if nargout > 14
    mumat = zeros(maxiter,k,d);
    pmat = zeros(maxiter,k);
    vmat = pmat;
    Lmat = zeros(maxiter,11);
end
pmat = zeros(maxiter,k);
Lprev = -inf;

for kind = 1:k
    
    lambdatmp = lambda(:,:,kind);
    
    %     for nind = 1:n
    %         mandist(nind,kind) = X(nind,:)*lambdatmp*X(nind,:)' - ...
    %             2*X(nind,:)*lambdatmp*mu(kind,:)' + ...
    %             trace(squeeze(mumu(:,:,kind))*lambdatmp);
    %     end
    %
    
    
    %
    %         for nind = 1:n
    %             mandist2(nind,kind) = X(nind,:)*lambdatmp*X(nind,:)' - ...
    %                 2*X(nind,:)*lambdatmp*mu(kind,:)' + ...
    %                 trace(mumu(:,:,kind)*lambdatmp);
    %
    %             a2(nind,kind) = (v(kind) + s(nind,kind)*d)/2;
    %             b2(nind,kind) = (v(kind) + s(nind,kind)*mandist2(nind,kind))/2;
    %         end
    
    mulambda = repmat(mu(kind,:)',1,d) .* lambdatmp;
    d1 = Xprod * lambdatmp(:);
    d2 = X2 * mulambda(:);
    
    mandist(:,kind) = d1 - 2*d2 + trace(mumu(:,:,kind)*lambdatmp);
    
    logdetlambda(kind) = log(det(lambda(:,:,kind)));
    
end
%
% s = gamma0;
% v = v0;
%
% for kind = 1:k
%
%     lambdatmp = lambda(:,:,kind);
%     %
%     %         for nind = 1:n
%     %             mandist2(nind,kind) = X(nind,:)*lambdatmp*X(nind,:)' - ...
%     %                 2*X(nind,:)*lambdatmp*mu(kind,:)' + ...
%     %                 trace(mumu(:,:,kind)*lambdatmp);
%     %
%     %             a2(nind,kind) = (v(kind) + s(nind,kind)*d)/2;
%     %             b2(nind,kind) = (v(kind) + s(nind,kind)*mandist2(nind,kind))/2;
%     %         end
%
%     lambda1 = shiftdim(repmat(lambdatmp,[1,1,n]),2);
%     mu1 = repmat(mu(kind,:),[n,1,d]);
%
%     d1 = sum(sum(Xprod .* lambda1,2),3);
%     d2 = sum(sum(mu1 .* X2 .* lambda1,2),3);
%
%     mandist(:,kind) = d1 - 2*d2 + trace(mumu(:,:,kind)*lambdatmp);
%
%     a(:,kind) = (v(kind) + s(:,kind)*d)/2;
%     b(:,kind) = (v(kind) + s(:,kind).*mandist(:,kind))/2;
%
%
% end
%
% u = a ./ b;
% logu = psi(a) - log(b);
%
u = rand(n,k)/2+1/2;
logu = log(u);

r = zeros(n,k);

onesvec = ones(n,1);

pidimconst = d*log(2*pi);

pprev = p0;
muprev = mu0;

iter = 0;
conv = 0;

logC_W = -eta0*.5*log(det(W0)) - (eta0*d*.5*log(2) + d*(d-1)*.25*log(pi) + sum(gammaln((eta0+1-[1:d])/2)));

R = zeros(d,d,k);

Lterms0 = zeros(11,1);

vprev = v0;
lambdaprev = lambda0;

Lvc = zeros(maxiter,1);
LLvc = Lvc;
KLmat = zeros(maxiter,4);


while conv == 0 && iter < maxiter
    iter = iter + 1;
    
    for kind = 1:k
        
        r(:,kind) = max(exp((.5*logdetlambda(kind)-.5*pidimconst)*onesvec + ...
            d*.5*logu(:,kind) - .5*u(:,kind).* mandist(:,kind)),eps);
        
        %         for nind = 1:n
        %             r(nind,kind)  = max(exp(logp(kind)+.5*logdetlambda(kind) + d*.5*logu(nind,kind) - ...
        %                 .5*u(nind,kind)*mandist(nind,kind) - .5*d*log(2*pi)),eps);
        %         end
    end
    
    [alphahat,betahat,c,s2sm] = baumwelch(A,r',Pi);
    s = (alphahat .* betahat ./ repmat(c,k,1))';
    s = max(s,eps);
    
    
    eta = eta0 + sum(s);
    
    %     w = s;
    
    w = s .* u;
    
    
    for kind = 1:k
        
        %         Winv(:,:,kind) = W0inv;
        %         for nind = 1:n
        %             Winv(:,:,kind) = Winv(:,:,kind) + s(nind,kind)*u(nind,kind)*(X(nind,:)'*X(nind,:) -...
        %                 X(nind,:)'*mu(kind,:) - mu(kind,:)'*X(nind,:) + mumu(:,:,kind));
        %         end
        
        wtmp = repmat(sqrt(w(:,kind)),1,d);
        Xtmp1 = X .* wtmp;
        Xtmp2 = Xtmp1 .* wtmp;
        matmp = Xtmp2'*repmat(mu(kind,:),n,1);
        %         mumutmp = squeeze(sum(shiftdim(repmat(mumu(:,:,kind),[1,1,n]),2) .*  repmat(w(:,kind),[1,d,d])));
        %
        
        mumutmp = mumu(:,:,kind)*sum(w(:,kind));
        Winv(:,:,kind) = W0inv + Xtmp1'*Xtmp1 + mumutmp - matmp - matmp';
        
        W(:,:,kind) = inv(Winv(:,:,kind));
        lambda(:,:,kind) = eta(kind)*W(:,:,kind);
        
        logdetlambda(kind) = log(det(lambda(:,:,kind)));
        
        shat = sum(s(:,kind));
        
        vconst = 1 + sum(s(:,kind).*(logu(:,kind)-u(:,kind)))/shat;
        %         for nind = 1:n
        %             vconst = vconst + s(nind,kind)*(logu(nind,kind)-u(nind,kind))/shat;
        %         end
        
        v(kind) = abs(fzero(@(v_k) vfun(v_k,vconst),vprev(kind)));
        
        R(:,:,kind) = lambda(:,:,kind) * sum(w(:,kind)) + p0*eye(d);
        mtmp = inv(R(:,:,kind)) * lambda(:,:,kind)*sum(repmat(w(:,kind),1,d).*X)';
        mu(kind,:) = mtmp;
        mumu(:,:,kind) = mtmp*mtmp' + inv(R(:,:,kind));
    end
    
    
    for kind = 1:k
        
        lambdatmp = lambda(:,:,kind);
        %
        %         for nind = 1:n
        %             mandist2(nind,kind) = X(nind,:)*lambdatmp*X(nind,:)' - ...
        %                 2*X(nind,:)*lambdatmp*mu(kind,:)' + ...
        %                 trace(mumu(:,:,kind)*lambdatmp);
        %
        %             a2(nind,kind) = (v(kind) + s(nind,kind)*d)/2;
        %             b2(nind,kind) = (v(kind) + s(nind,kind)*mandist2(nind,kind))/2;
        %         end
        
        %         lambda1 = shiftdim(repmat(lambdatmp,[1,1,n]),2);
        %         mu1 = repmat(mu(kind,:),[n,1,d]);
        %
        %         d1 = sum(sum(Xprod .* lambda1,2),3);
        
        mulambda = repmat(mu(kind,:)',1,d) .* lambdatmp;
        d1 = Xprod * lambdatmp(:);
        d2 = X2 * mulambda(:);
        
        %         d2 = sum(sum(mu1 .* X2 .* lambda1,2),3);
        %
        mandist(:,kind) = d1 - 2*d2 + trace(mumu(:,:,kind)*lambdatmp);
        
        a(:,kind) = (v(kind) + s(:,kind)*d)/2;
        b(:,kind) = (v(kind) + s(:,kind).*mandist(:,kind))/2;
        
        
    end
    
    
    u = a ./ b;
    
    psia = psi(a);
    logb = log(b);
    
    logu = psia - logb;
    
    Aprev = A;
    
    alpha = s2sm + alphak;
    alphasm = sum(s(1:n-1,:)) + alpha0;
    alphasmat = repmat(alphasm',1,k);
    
    logA = psi(alpha) - psi(alphasmat);
    A = exp(logA);
    
    Pitmp = s(1,:)' + alphak;
    Pismtmp = sum(s(1,:)') + alpha0;
    
    logPi = psi(Pitmp) - psi(Pismtmp);
    
    Pi = exp(logPi);
    
    %     KL_A = -sum(gammaln(alphasm) - gammaln(alpha0) - sum(gammaln(alpha) - gammaln(alphak) - (s2sm .* logA)));
    %     KL_Pi = -(gammaln(Pismtmp) - gammaln(alpha0) - sum(gammaln(Pitmp) - gammaln(alphak) - Pitmp .* logPi));
    %
    %     KL_mu = 0;
    %     KL_lambda = 0;
    %
    %     for kind = 1:k
    %         Rk = squeeze(R(:,:,kind));
    %         s_smk = s_sm(kind);
    %         logq_muk = -.5*(log(det(Rk)));
    %         logp_muk = -.5*(d*log(p0) - p0*norm(mu(kind,:)));
    %
    %         KL_mu = KL_mu + logp_muk - logq_muk;
    %
    %
    %         etatmp = eta(kind);
    %         Wtmp = W(:,:,kind);
    %
    %         logC_Wm = -etatmp*.5*log(det(Wtmp)) - (etatmp*d*.5*log(2) + d*(d-1)*.25*log(pi) + sum(gammaln((etatmp+1-[1:d])/2)));
    %
    %         logq_lambdak = logC_Wm + .5*(eta(kind) - d - 1)*logdetlambda(kind) - .5*eta(kind)*d;
    %         logp_lambdak = logC_W + .5*(eta0 - d - 1)*logdetlambda(kind) - .5*trace(W0inv*lambda(:,:,kind));
    %
    %         KL_lambda = KL_lambda + logp_lambdak - logq_lambdak;
    %
    %
    %     end
    %
    %     KL = KL_A + KL_mu + KL_lambda + KL_Pi;
    LL = -sum(log(c));
    
    %     Lvc(iter) = L;
    %     KLmat(iter,:) = [KL_A KL_mu KL_lambda KL_Pi];
    LLvc(iter) = LL;
    
    if nargout > 15
        mumat(iter,:,:) = mu;
        pmat(iter,:) =  sum(s) / sum(s(:));
        vmat(iter,:) = v;
    end
    
    
    pmat(iter,:) =  sum(s) / sum(s(:));
    Adiff = (A(:) - Aprev(:));
    mdiff = (mu - muprev);
    lambdadiff = lambda - lambdaprev;
    vdiff = (v - vprev) ./ vprev;
    %     Ldiff = 100*(L - Lprev) / Lprev;
    
    if ~Lconvflg
        conv = max(abs(Adiff)) < 0.001 && max(abs(mdiff(:))) < 0.01 && max(abs(lambdadiff(:))) < 0.001 && max(abs(vdiff)) < 0.1;
    else
        conv = abs(Ldiff) < 0.01;
    end
    
    Aprev = A;
    muprev = mu;
    lambdaprev = lambda;
    vprev = v;
    %
    %     Lprev = L;
    
    
    
end

Lvc = Lvc(1:iter);
LLvc = LLvc(1:iter);
KLmat = KLmat(1:iter,:);
pmat = pmat(1:iter,:);

sigma = zeros(d,d,k);

DIC = 0;

for kind = 1:k
    sigma(:,:,kind) = inv(lambda(:,:,kind));
end


if nargout > 14
    Lmat = Lmat(1:iter,:);
    mumat = mumat(1:iter,:,:);
    pmat = pmat(1:iter,:);
    vmat = vmat(1:iter,:);
end

L = LL;

%
% term_A = sum(sum(s2sm .* logA));
% term_mu = 0;
% term_lambda = 0;
% term_Pi = sum(Pitmp .* log(max(Pi,eps)));
%
% for kind = 1:k
%     Rk = squeeze(R(:,:,kind));
%     s_smk = sum(s(:,kind));
%     logq_muk = .5*(log(det(Rk)));
%     logp_muk = -.5*(d*log(p0) - p0*norm(mu(kind,:)));
%
%     term_mu = term_mu + s_smk*(logq_muk - logp_muk);
%
%
%     etatmp = eta(kind);
%     Wtmp = W(:,:,kind);
%
%     logC_Wm = -etatmp*.5*log(det(Wtmp)) - (etatmp*d*.5*log(2) + d*(d-1)*.25*log(pi) + sum(gammaln((etatmp+1-[1:d])/2)));
%
%     logq_lambdak = logC_Wm + .5*(eta(kind) - d - 1)*logdetlambda(kind) - .5*eta(kind)*d;
%     logp_lambdak = logC_W + .5*(eta0 - d - 1)*logdetlambda(kind) - .5*trace(W0inv*lambda(:,:,kind));
%
%     term_lambda = term_lambda + s_smk*(logq_lambdak - logp_lambdak);
%
%
% end
%
% pD = -2*(term_A + term_mu + term_lambda + term_Pi);
%
% term2_A = sum(sum(s2sm .* log(max(alpha ./ (alphasmat'),eps))));
% % term2_mu =


%
%     for kind = 1:k
%
%         Lterms(1) = Lterms(1) + ...
%             .5*sum(s(:,kind) .* (logdetlambda(kind) - pidimconst + d*logu(:,kind) - u(:,kind) .* mandist(:,kind)));
%
%         Lterms(2) = Lterms(2) + .5*(d*log(p0/(2*pi)) - p0*norm(mu(kind,:)));
%
%         Lterms(3)= Lterms(3) + ...
%             logC_W + .5*(eta0 - d - 1)*logdetlambda(kind) - .5*trace(W0inv*lambda(:,:,kind));
%
%         Lterms(4) = Lterms(4) + n*(.5*v(kind)*log(.5*v(kind)) - gammaln(.5*v(kind))) + ...
%             sum((.5*v(kind)-1)*logu(:,kind) - .5*v(kind)*u(:,kind));
%
%         Lterms(5) = Lterms(5) + (alphak-1)*sum(logA(kind,:)) - gammaln(alphak);
%
%         Lterms(6) = Lterms(6) + sum(s(:,kind)*sum(logA(kind,:)));
%
%         Lterms(7) = Lterms(7) -(.5*log(det(R(:,:,kind))) - .5*d*(1+log(2*pi)));
%
%         etatmp = eta(kind);
%         Wtmp = W(:,:,kind);
%
%         logC_Wm = -etatmp*.5*log(det(Wtmp)) - (etatmp*d*.5*log(2) + d*(d-1)*.25*log(pi) + sum(gammaln((etatmp+1-[1:d])/2)));
%
%         Lterms(8) = Lterms(8) - (logC_Wm + .5*(eta(kind) - d - 1)*logdetlambda(kind) - .5*eta(kind)*d);
%
%         Lterms(9) = Lterms(9) - ...
%             sum((a(:,kind)-1).*psia(:,kind) - a(:,kind) - gammaln(a(:,kind)) + logb(:,kind));
%
%         Lterms(10) = Lterms(10) - (sum((alpha(:,kind)-1)*logA(kind,:)) - gammaln(sum(alpha(:,kind))));
%
%         Lterms(11) = Lterms(11) - sum(s(:,kind) .* log(s(:,kind)));
%     end
%
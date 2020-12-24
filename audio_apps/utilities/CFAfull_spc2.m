function [W,psi,phi,sigma,omega,lambda,iter,n_fail,logp,Chat] = CFAfull_spc2(X,Q,bootN,m,D,R,V)

if nargin < 4
    m = 1;
end

[n,d] = size(X);
Qd = size(Q,2);

if nargin < 5
    D = inv(tril(ones(d)));
end

if nargin < 6
    R = D*Q;
end

if nargin < 7
    V = zeros(d,0);
end

Dd = size(D,2);
Rd = size(R,2);
Vd = size(V,2);

W = zeros(d,m,bootN);
sigma = zeros(Dd,bootN);
psi = zeros(d,bootN);
phi = zeros(Qd,bootN);
omega = zeros(Rd,bootN);
lambda = zeros(Vd,bootN);
iter = zeros(1,bootN);
logp = iter;

W0 = W;
psi0 = psi;
sigma0 = sigma;
phi0 = phi;
omega0 = omega;
lambda0 = lambda;


C = cov(X)*n/(n-1);
vr = diag(C);
vrmn = mean(vr);

i = 0;
n_fail = 0;

% h = waitbar(0/bootN,'Calculating FA parameters');

while i < bootN

    %     indtmp = ceil(n*rand(n,1));
    
    
    W0tmp = repmat(sqrt(vr),1,m).*rand(d,m);
    
    psi0tmp = vr.*rand(d,1);
    
    sigma0tmp = .5*rand(Dd,1)*vrmn;
    phi0tmp = rand(Qd,1)*vrmn;
    omega0tmp = .5*rand(Rd,1)*vrmn;
    lambda0tmp = .5*rand(Vd,1)*vrmn;
    
%     sigma0tmp = diag(diag(D'*diag(vr.*rand(d,1))*D ./ diag(diag(D'*D))));
%     phi0tmp = diag(diag(Q'*diag(vr.*rand(d,1))*Q ./ diag(diag(Q'*Q))));
%     omega0tmp = diag(diag(R'*diag(vr.*rand(d,1))*R ./ diag(diag(R'*R))));
    
    [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,lambdatmp,itertmp,degflg] = CFA_full2(X,Q,m,W0tmp,sigma0tmp,psi0tmp,phi0tmp,omega0tmp,lambda0tmp,D,R,V,0.01);

    if itertmp < 10000 && degflg
        i = i + 1;

        psi(:,i) = psitmp;
        psi0(:,i) = psi0tmp;
    
        if m > 0
            W(:,:,i) = Wtmp;
            W0(:,:,i) = W0tmp;           
        end
        
        if Dd > 0
            sigma(:,i) = sigmatmp;
            sigma0(:,i) = sigma0tmp;
        end
        
        if Qd > 0
            phi(:,i) = phitmp;
            phi0(:,i) = phi0tmp;
        end
        
        if Rd > 0
            omega(:,i) = omegatmp;
            omega0(:,i) = omega0tmp;
        end
        
        if Vd > 0
            lambda(:,i) = lambdatmp;
            lambda0(:,i) = lambda0tmp;
        end

        iter(i) = itertmp;
        
        Chat = Wtmp*Wtmp' + D*diag(sigmatmp)*D' + Q*diag(phitmp)*Q' + diag(psitmp) + R*diag(omegatmp)*R' + V*diag(lambdatmp)*V';
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
%         h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

% [logp,maxind] = max(logp);
% 
% sigma = sigma(:,maxind);
% W = W(:,:,maxind);
% psi = psi(:,maxind);
% omega = omega(:,maxind);
% lambda = lambda(:,maxind);
% phi = phi(:,maxind);
% iter = iter(maxind);
% 
% 
% Chat = W*W' + D*diag(sigma)*D' + Q*diag(phi)*Q' + diag(psi) + R*diag(omega)*R' + V*diag(lambdatmp)*V';

inds = find(logp >= prctile(logp,95));

logp = logp(inds);
W0 = W0(:,:,inds);
psi0 = psi0(:,inds);
phi0 = phi0(:,inds);
sigma0 = sigma0(:,inds);
omega0 = omega0(:,inds);
lambda0 = lambda0(:,inds);

logprev = 10*min(logp);

for i = 1:length(inds)
    [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,lambdatmp,itertmp,degflg] = CFA_full2(X,Q,m,W0(:,:,i),sigma0(:,i),psi0(:,i),phi0(:,i),omega0(:,i),lambda0(:,i),D,R,V,0.0001);
    
    Chatmp = Wtmp*Wtmp' + D*diag(sigmatmp)*D' + Q*diag(phitmp)*Q' + diag(psitmp) + R*diag(omegatmp)*R' + V*diag(lambdatmp)*V';
    logptmp = -n/2*(d*log(2*pi) + log(det(Chatmp)) + trace(inv(Chatmp)*C));
    
    if logptmp > logprev
        W = Wtmp;
        psi = psitmp;
        sigma = sigmatmp;
        phi = phitmp;
        omega = omegatmp;
        lambda = lambdatmp;
        Chat = Chatmp;
        logp = logptmp;
        iter = itertmp;
        
    end
    
end
% 
% close h
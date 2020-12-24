function [W,psi,phi,sigma,omega,iter,n_fail,logp,Chat,W0,psi0,phi0,sigma0,omega0] = CFAfull_spc_waitbar(X,Q,initN,m,D,R)

if nargin < 4
    m = 1;
end

[n,d] = size(X);
Qd = size(Q,2);

if nargin < 5
    D = -diff(eye(d))';
end

if nargin < 6
    R = Q;
end

Dd = size(D,2);
Rd = size(R,2);

W = zeros(d,m,initN);
sigma = zeros(Dd,initN);
psi = zeros(d,initN);
phi = zeros(Qd,initN);
omega = zeros(Rd,initN);
iter = zeros(1,initN);
logp = iter;

W0 = W;
psi0 = psi;
sigma0 = sigma;
phi0 = phi;
omega0 = omega;


C = cov(X)*(n-1)/n;
vr = diag(C);
vrmn = mean(vr);

i = 0;
n_fail = 0;

h = waitbar(0/initN,'Calculating FA parameters');

while i < initN
    
    %     indtmp = ceil(n*rand(n,1));
    
    
    W0tmp = repmat(sqrt(vr),1,m).*rand(d,m);
    
    psi0tmp = vr.*rand(d,1);
    
    sigma0tmp = .5*rand(Dd,1)*vrmn;
    phi0tmp = rand(Qd,1)*vrmn;
    omega0tmp = .5*rand(Rd,1)*vrmn;
    
    %     sigma0tmp = diag(diag(D'*diag(vr.*rand(d,1))*D ./ diag(diag(D'*D))));
    %     phi0tmp = diag(diag(Q'*diag(vr.*rand(d,1))*Q ./ diag(diag(Q'*Q))));
    %     omega0tmp = diag(diag(R'*diag(vr.*rand(d,1))*R ./ diag(diag(R'*R))));
    
    %     [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,itertmp,degflg] = CFA_full(X,Q,m,W0tmp,sigma0tmp,psi0tmp,phi0tmp,omega0tmp,D,R,.01,1);
    [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,itertmp,degflg] = CFA_full(X,Q,m,W0tmp,sigma0tmp,psi0tmp,phi0tmp,omega0tmp,D,R,.1,1);

    if itertmp < 100000 && degflg
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
        
        iter(i) = itertmp;
        
        Chat = Wtmp*Wtmp' + D*diag(sigmatmp)*D' + Q*diag(phitmp)*Q' + diag(psitmp) + R*diag(omegatmp)*R';
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));

        
        h = waitbar(i/initN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end
    
end

close(h)

% [logp,maxind] = max(logp);
%
% sigma = sigma(:,maxind);
% W = W(:,:,maxind);
% psi = psi(:,maxind);
% omega = omega(:,maxind);
% phi = phi(:,maxind);
% iter = iter(maxind);
% Chat = W*W' + D*diag(sigma)*D' + Q*diag(phi)*Q' + diag(psi) + R*diag(omega)*R';

inds = find(logp >= prctile(logp,95));
% inds = find(logp == max(logp));

logp = logp(inds);
W0all = W0(:,:,inds);
psi0all = psi0(:,inds);
phi0all = phi0(:,inds);
sigma0all = sigma0(:,inds);
omega0all = omega0(:,inds);

logprev = 10*min(logp);

for i = 1:length(inds)
    
    W0tmp = W0all(:,:,i);
    sigma0tmp = sigma0all(:,i);
    psi0tmp = psi0all(:,i);
    phi0tmp = phi0all(:,i);
    omega0tmp = omega0all(:,i);
    
    [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,itertmp,degflg,z,u] = CFA_full(X,Q,m,W0tmp,sigma0tmp,psi0tmp,phi0tmp,omega0tmp,D,R,0.1);
    
    Chatmp = Wtmp*Wtmp' + D*diag(sigmatmp)*D' + Q*diag(phitmp)*Q' + diag(psitmp) + R*diag(omegatmp)*R';
    logptmp = -n/2*(d*log(2*pi) + log(det(Chatmp)) + trace(inv(Chatmp)*C));

    if logptmp > logprev
        W = Wtmp;
        psi = psitmp;
        sigma = sigmatmp;
        phi = phitmp;
        omega = omegatmp;
        Chat = Chatmp;
        logp = logptmp;
        iter = itertmp;
        
        logprev = logptmp;
        
        W0 = W0tmp;
        psi0 = psi0tmp;
        sigma0 = sigma0tmp;
        phi0 = phi0tmp;
        omega0 = omega0tmp;
        
    end
    
end


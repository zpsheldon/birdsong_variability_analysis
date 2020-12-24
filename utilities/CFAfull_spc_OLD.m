function [W,psi,phi,sigma,omega,iter,n_fail,logp,Chat] = CFAfull_spc_OLD(X,Q,bootN,m,D,R)

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

Dd = size(D,2);
Rd = size(R,2);

W = zeros(d,m,bootN);
sigma = zeros(Dd,bootN);
psi = zeros(d,bootN);
phi = zeros(Qd,bootN);
omega = zeros(Rd,bootN);
iter = zeros(1,bootN);
logp = iter;

W0 = W;
psi0 = psi;
sigma0 = sigma;
phi0 = phi;
omega0 = omega;




C = cov(X)*n/(n-1);
vr = diag(C);
vrmn = mean(vr);

i = 0;
n_fail = 0;

%h = waitbar(0/bootN,'Calculating FA parameters');

while i < bootN

    %     indtmp = ceil(n*rand(n,1));
    
    
    W0tmp = repmat(sqrt(vr),1,m).*rand(d,m);
    
    psi0tmp = diag(vr.*rand(d,1));
    
    sigma0tmp = diag(2*rand(Dd,1)*vrmn);
    phi0tmp = diag(rand(Qd,1)*vrmn);
    omega0tmp = diag(2*rand(Rd,1)*vrmn);
%     
%     sigma0tmp = diag(diag(D'*diag(vr.*rand(d,1))*D ./ diag(diag(D'*D))));
%     phi0tmp = diag(diag(Q'*diag(vr.*rand(d,1))*Q ./ diag(diag(Q'*Q))));
%     omega0tmp = diag(diag(R'*diag(vr.*rand(d,1))*R ./ diag(diag(R'*R))));
    
    [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,itertmp,degflg] = CFA_full(X,Q,m,W0tmp,sigma0tmp,psi0tmp,phi0tmp,omega0tmp,D,R);

    if itertmp < 10000 && degflg
        i = i + 1;

        W(:,:,i) = Wtmp;
        sigma(:,i) = diag(sigmatmp);
        psi(:,i) = diag(psitmp);
        
        W0(:,:,i) = W0tmp;
        sigma0(:,i) = diag(sigma0tmp);
        psi0(:,i) = diag(psi0tmp);
    
        if Qd > 0
            phi(:,i) = diag(phitmp);
            phi0(:,i) = diag(phi0tmp);
            
            omega(:,i) = diag(omegatmp);
            omega0(:,i) = diag(omega0tmp);
        end

        iter(i) = itertmp;
        
        Chat = Wtmp*Wtmp' + D*sigmatmp*D' + Q*phitmp*Q' + psitmp + R*omegatmp*R';
        logp(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
        
%         h = waitbar(i/bootN,h,'Calculating FA parameters');
        
    else
        n_fail = n_fail + 1
    end

end

[logp,maxind] = max(logp);
W = W(:,:,maxind);
psi = psi(:,maxind);
phi = phi(:,maxind);
sigma = sigma(:,maxind);
omega = omega(:,maxind);

Chat = W*W' + D*diag(sigma)*D' + Q*diag(phi)*Q' + diag(psi) + R*diag(omega)*R';

%close(h)

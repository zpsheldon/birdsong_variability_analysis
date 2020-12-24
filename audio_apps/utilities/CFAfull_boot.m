function [bootstrct,Wall,psiall,phiall,sigmall,omegall,logpall] = CFAfull_boot(X,Q,bootN,D,R,What,psihat,phihat,sigmahat,omegahat)

[n,d] = size(X);
Qd = size(Q,2);
Dd = size(D,2);
Rd = size(R,2);

m = size(What,2);

C = cov(X)*(n-1)/n;
vr = diag(C);
vrmn = mean(vr);

i = 0;
n_fail = 0;

h = waitbar(0/bootN,'Bootstrapping FA parameters');

Wall = zeros(d,m,bootN);
sigmall = zeros(Dd,bootN);
psiall = zeros(d,bootN);
phiall = zeros(Qd,bootN);
omegall = zeros(Rd,bootN);
logpall = zeros(1,bootN);

% P = zeros(n,bootN);

Chat = What*What' + D*diag(sigmahat)*D' + Q*diag(phihat)*Q' + diag(psihat) + R*diag(omegahat)*R';

% X_e = X * chol(Chat) * inv(chol(C));

X_e = X;

i = 0;

while i < bootN
   
    bootinds = ceil(rand(1,n)*n);
    Xtmp = X_e(bootinds,:);
   
%     [Wtmp,sigmatmp,psitmp,phitmp,omegatmp,itertmp,degflg] = CFA_full(Xtmp,Q,m,W0,sigma0,psi0,phi0,omega0,D,R,0.0001);
%     
    [Wtmp,psitmp,phitmp,sigmatmp,omegatmp,itertmp] = CFAfull_spc(Xtmp,Q,100,m,D,R,0);
   
    if itertmp < 100000 %&& degflg
        i = i + 1;
               
%         tabtmp = tabulate(bootinds);
%         P(tabtmp(:,1),i) = tabtmp(:,3);

        psiall(:,i) = psitmp;
       
        if m > 0
            Wall(:,:,i) = rotate_maxaccum(Wtmp,What);
        end
       
        if Dd > 0
            sigmall(:,i) = sigmatmp;
        end
       
        if Qd > 0
            phiall(:,i) = phitmp;
        end
       
        if Rd > 0
            omegall(:,i) = omegatmp;
        end

        Chat = Wtmp*Wtmp' + D*diag(sigmatmp)*D' + Q*diag(phitmp)*Q' + diag(psitmp) + R*diag(omegatmp)*R';
        logpall(i) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*C));
       
        h = waitbar(i/bootN,h,'Bootstrapping FA parameters');
       
    else
        n_fail = n_fail + 1
    end
   
   
end


bootstrct.sqrtpsimn = mean(sqrt(psiall),2);
bootstrct.sqrtpsisem = std(sqrt(psiall),0,2);
bootstrct.psiall = psiall;

bootstrct.Wmn = [];
bootstrct.Wsem = [];
bootstrct.sqrtphimn = [];
bootstrct.sqrtphisem = [];
bootstrct.sqrtomegamn = [];
bootstrct.sqrtomegasem = [];
bootstrct.sqrtsigmamn = [];
bootstrct.sqrtsigmasem = [];

bootstrct.Wall = [];
bootstrct.phiall = [];
bootstrct.omegall = [];
bootstrct.sigmall = [];

if m > 0
    bootstrct.Wmn = mean(Wall,3);
    bootstrct.Wsem = std(Wall,0,3);
    bootstrct.Wall = Wall;
end

if Qd > 0
    bootstrct.sqrtphimn = mean(sqrt(phiall),2);
    bootstrct.sqrtphisem = std(sqrt(phiall),0,2);
    bootstrct.phiall = phiall;
end

if Rd > 0
    bootstrct.sqrtomegamn = mean(sqrt(omegall),2);
    bootstrct.sqrtomegasem = std(sqrt(omegall),0,2);
    bootstrct.omegall = omegall;
end

if Dd > 0
    bootstrct.sqrtsigmamn = mean(sqrt(sigmall),2);
    bootstrct.sqrtsigmasem = std(sqrt(sigmall),0,2);
    bootstrct.sigmall = sigmall;
end

bootstrct.logpmn = mean(logpall);
bootstrct.logpsem = std(logpall);
bootstrct.logpall = logpall;

close(h)
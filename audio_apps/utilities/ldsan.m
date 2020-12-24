function [A,C,Gamma,Sigma,expz,expz2sm,Avc] = ldsan(X,dt)

[n,d] = size(X);

if nargin < 2
    dt = ones(1,n-1);
end

m = 1 + 2*d-1;
D = -diff(eye(d))';

P0 = rand(m,m);
mu0 = rand(m,1);
A = diag(diag(rand(m)));
% A = diag(.5*ones(m,1));
% A = [.5 0;0 .5];
C = [rand(d,1) eye(d) D];

Gamma = [1 zeros(1,2*d-1);zeros(2*d-1,1) diag(rand(2*d-1,1))];
% Sigma = Gamma(2:2+d-1,2:2+d-1);

Sigma = diag(rand(d,1));

mu = zeros(n,m);
V = zeros(n,m,m);
P = zeros(n,m,m);

if m > 1
    K = zeros(n,m,d);
else
    K = zeros(n,d);
end

J = V;
muhat = mu;
Vhat = V;

maxiter = 20;
convflg = 0;
iter = 0;

Gammatmp = zeros(m);

while iter <= maxiter && ~convflg
    
    iter = iter + 1
    
    K(1,:,:) = P0 * C' * inv(C*P0*C' + Sigma);
    
    mu(1,:) = mu0 + squeeze(K(1,:,:)) * (X(1,:)' - C*mu0);
    V(1,:,:) = (eye(m) - squeeze(K(1,:,:))*C)*P0;
    
    Gammatmp(:) = 0;
    
    for i = 2:n
        Atmp = A ^ dt(i-1);
        kvc = 2*[0:dt(i-1)-1];
        
        Gammatmp(:) = 0;
        for j = 1:length(kvc)
            Gammatmp = Gammatmp + Gamma * (A ^ kvc(j));
        end
        
        P(i-1,:,:) = Atmp * squeeze(V(i-1,:,:)) * Atmp' + Gammatmp;
        K(i,:,:) = squeeze(P(i-1,:,:)) * C' * inv(C * squeeze(P(i-1,:,:)) * C' + Sigma);
        mu(i,:) = Atmp*mu(i-1,:)' + squeeze(K(i,:,:)) * (X(i,:)' - C * Atmp * mu(i-1,:)');
        V(i,:,:) = (eye(m) - squeeze(K(i,:,:))*C) * squeeze(P(i-1,:,:));
    end
    
    muhat(n,:) = mu(n,:);
    Vhat(n,:,:) = V(n,:,:);
    
    for i = n-1:-1:1
        Atmp = A ^ dt(i);
        J(i,:,:) = squeeze(V(i,:,:)) * Atmp' * inv(squeeze(P(i,:,:)));
        muhat(i,:) = mu(i,:)' + squeeze(J(i,:,:)) * (muhat(i+1,:)' - Atmp * mu(i,:)');
        Vhat(i,:,:) = squeeze(V(i,:,:)) + squeeze(J(i,:,:)) * (squeeze(Vhat(i+1,:,:)) - squeeze(P(i,:,:))) * squeeze(J(i,:,:))';
    end
    
    expz2 = zeros(n,m,m);
    expzzmin1 = expz2;
%     expzzmine = expz2;
    
    expz = muhat;
    
    for i = 2:n
        expzzmin1(i,:,:) = squeeze(Vhat(i,:,:))*squeeze(J(i-1,:,:))' + muhat(i,:)'*muhat(i-1,:);
%         expzzmine(i,:,:) = squeeze(expzzmin1(i,:,:))^(1/dt(i-1));
    end
    
    for i = 1:n
        expz2(i,:,:) = squeeze(Vhat(i,:,:)) + muhat(i,:)'*muhat(i,:);
    end
    
    mu0 = expz(1,:)';
    P0 = squeeze(expz2(1,:,:)) - expz(1,:)'*expz(1,:);
    
    Aold = A;
%     A = diag(diag(squeeze(sum(expzzmin1(2:end,:,:))) * inv(squeeze(sum(expz2(1:end-1,:,:))))));

    A = squeeze(sum(expzzmin1(2:end,:,:))) * inv(squeeze(sum(expz2(1:end-1,:,:))));
    A = diag(diag(A));
    
%     invexpzzmin1sm = inv(squeeze(sum(expzzmin1(2:end,:,:))));
%     invexpz2sm = inv(squeeze(sum(expz2(1:end-1,:,:))));
%     coeffs = zeros(n-1,m,m);
%     
%     for i = 1:n-1
%         coeffs(i,:,:) = Aold^(dt(i)-1)*squeeze(expzzmin1(i+1,:,:))* invexpz2sm;
%     end
%     
%     A = diag(diag(squeeze(sum(coeffs))));
    
%     
%     ftmp = @(y,k,c1,c2)((y.^k)*c1-c2);
% 
%     
%     for j = 1:m
%         A(j,j) = fzero(@(x) ftmp(x,dt,squeeze(coeffs(:,j,j)),1),Aold(j,j));
%     end
    
    Avc(iter)= A(1,1);
    A(1,1)
    
    if max(abs(diag(A-Aold)))<0.0001
        convflg = 1;
    end
    
    Gammatmp = zeros(size(Gamma));
    
    for i = 1:n-1
        Atmp = A ^ dt(i);
        Gammatmp = Gammatmp + squeeze(expz2(i+1,:,:)) + Atmp*squeeze(expzzmin1(i+1,:,:)) + squeeze(expzzmin1(i+1,:,:))*Atmp + Atmp * squeeze(expz2(i,:,:)) * Atmp';
    end
    Gammatmp = Gammatmp / (n-1);
    
%     Gamma = [1 zeros(1,2*d-1);zeros(2*d-1,1) diag(diag(Gammatmp(2:end,2:end)))];
    
    C = X'*expz * inv(squeeze(sum(expz2)));
    C = [C(:,1) eye(d) D];
    
    Sigma = zeros(d);
    
    for i = 1:n
        Sigma = Sigma + X(i,:)'*X(i,:) - C*expz(i,:)'*X(i,:) - X(i,:)'*expz(i,:)*C' + C*squeeze(expz2(i,:,:))*C';
    end
    
    Sigma = Sigma / n;
    
    %     xhat = expz*C';
    %
    %     Sigma = covx - 2 * xhat'*X/n + C*squeeze(sum(expz2))*C'/n;
    
%     Sigma = diag(diag(Sigma));
    
end

expz2sm = squeeze(sum(expz2));

function paramstrct = CFA_fullmod_noboot(X,seqarr,featarr,norptflg,m_glob,mflg)

if nargin<4
    norptflg = 0;
end

Q = mkQmat(seqarr);
R = mkRmat(featarr);
D = -diff(eye(length(seqarr)))';

Dnull = zeros(length(seqarr),0);

if nargin >= 5
    if nargin >= 6 && mflg
        mvc = m_glob;
    else
        mvc = 1:m_glob;
    end
else
    mvc = 1:3;
end

Darr = {D};
if ~size(Q,2) || norptflg
    Qarr = {Dnull};
else
    Qarr = {Dnull,Q};
end

if ~size(R,2) || norptflg
    Rarr = {Dnull};
else
    Rarr = {Dnull,R};
end

R = Dnull;
Rarr = {Dnull};

N = size(X,1);
C = cov(X);
Chat0 = diag(diag(C));
logp0 = -N/2*(size(X,2)*log(2*pi) + log(det(Chat0)) + trace(inv(Chat0)*C));

tag = {'-','+'};

% modelnm = length(mvc)*length(Darr)*length(Qarr)*length(Rarr);
modelnm = length(mvc)*length(Qarr)*length(Rarr);

paramstrct.iter = zeros(1,modelnm);
paramstrct.n_fail = paramstrct.iter;
paramstrct.k = paramstrct.iter;
paramstrct.logp = paramstrct.iter;

paramstrct.modelab{1} = 'no factors';
paramstrct.k(1) = size(X,2);
paramstrct.logp(1) = logp0;

modelind = 0;

for mind = 1:length(mvc)
    
    mtmp = mvc(mind);
    
    for Dind = 1:length(Darr)
        
        Dtmp = Darr{Dind};
        
        for Qind = 1:length(Qarr)
            
            Qtmp = Qarr{Qind};
            
            for Rind = 1:length(Rarr)
                
                Rtmp = Rarr{Rind};
                
                modelind = modelind + 1
                
                if modelind > 0
                    
                    [Wtmp,psitmp,phitmp,sigmatmp,omegatmp,itertmp,n_fail_tmp,logptmp,Shatmp,W0tmp,psi0tmp,phi0tmp,sigma0tmp,omega0tmp] = ...
                        CFAfull_spc(X,Qtmp,100,mtmp,Dtmp,Rtmp);
                    
                    ktmp = numel(Wtmp)+numel(psitmp)+numel(phitmp)+numel(sigmatmp)+numel(omegatmp);
                    
                    labtmp = [num2str(mtmp) ' glob, ' tag{Dind} 'jitter, ' tag{Qind} 'rpt. ns., ' tag{Rind} 'rpt. jttr.'];
                    
                    if mtmp > 0
                        paramstrct.W{modelind} = rotate_maxaccum(Wtmp);
                    end
                    
                    paramstrct.psi{modelind} = psitmp;
                    paramstrct.phi{modelind} = phitmp;
                    paramstrct.sigma{modelind} = sigmatmp;
                    paramstrct.omega{modelind} = omegatmp;
                    paramstrct.iter(modelind) = itertmp;
                    paramstrct.n_fail(modelind) = n_fail_tmp;
                    paramstrct.k(modelind) = ktmp;
                    paramstrct.logp(modelind) = logptmp;
                    paramstrct.modelab{modelind} = labtmp;
                    paramstrct.Shat{modelind} = Shatmp;
                    
                    paramstrct.W0{modelind} = W0tmp;
                    paramstrct.psi0{modelind} = psi0tmp;
                    paramstrct.phi0{modelind} = phi0tmp;
                    paramstrct.sigma0{modelind} = sigma0tmp;
                    paramstrct.omega0{modelind} = omega0tmp;
                    
                    paramstrct.Q{modelind} = Qtmp;
                    paramstrct.D{modelind} = Dtmp;
                    paramstrct.R{modelind} = Rtmp;
                    
                end
                
            end
            
        end
        
    end
    
end

paramstrct.bic = -2 * paramstrct.logp + log(N) * paramstrct.k;
paramstrct.aic = 2 * paramstrct.k - 2 * paramstrct.logp;

[minbic,modelind] = min(paramstrct.bic);

paramstrct.glob = paramstrct.W{modelind};
paramstrct.ind = sqrt(paramstrct.psi{modelind});
paramstrct.indrpt = sqrt(paramstrct.phi{modelind});
paramstrct.jitter = sqrt(paramstrct.sigma{modelind});
paramstrct.jitterpt = sqrt(paramstrct.omega{modelind});

paramstrct.iter = paramstrct.iter(modelind);
paramstrct.n_fail = paramstrct.n_fail(modelind);
paramstrct.k = paramstrct.k(modelind);
paramstrct.logp = paramstrct.logp(modelind);
paramstrct.modelab = paramstrct.modelab{modelind};
paramstrct.Shat = paramstrct.Shat{modelind};

paramstrct.Q = paramstrct.Q{modelind};
paramstrct.D = paramstrct.D{modelind};
paramstrct.R = paramstrct.R{modelind};

R = corrcoef(X);
Rhat = paramstrct.Shat ./ (std(X)'*std(X));

paramstrct.SRMR = sqrt(mean((vech(R)-vech(Rhat)).^2));


%
%
% [C_theta,semW,semN,semJ,semWvc,semNvc,semJvc] = asympcov_CFA2(paramstrct.W,diag(paramstrct.sqrtpsi.^2),diag(paramstrct.sqrtsigma.^2),X);

% paramstrct.Wsemhat = semWvc;
% paramstrct.sqrtpsisemhat = semNvc;
% paramstrct.sqrtsigmasemhat = semJvc;
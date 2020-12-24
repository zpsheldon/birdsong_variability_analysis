function [paramstrct_arr,permstrct] = CFA_fullperm(Xarr,seqarr,featarr,norptflg,m_glob,bootN)

if nargin<4
    norptflg = 0;
end

d = size(Xarr{1},2);

Xnm = length(Xarr);
Nsub = zeros(1,Xnm);
msub = Nsub;
Qdsub = Nsub;
Ddsub = Nsub;
Rdsub = Nsub;

for Xind = 1:Xnm
    Nsub(Xind) = size(Xarr{Xind},1);
    Carr{Xind} = cov(Xarr{Xind});
    paramstrct_arr{Xind} = CFA_fullmod_noboot(Xarr{Xind},seqarr,featarr,norptflg,m_glob);
    
    msub(Xind) = size(paramstrct_arr{Xind}.glob,2);
    
    Qdsub(Xind) = size(paramstrct_arr{Xind}.indrpt,1);
    Ddsub(Xind) = size(paramstrct_arr{Xind}.jitter,1);
    Rdsub(Xind) = size(paramstrct_arr{Xind}.jitterpt,1);
    
    Wall_arr{Xind} = zeros(d,msub(Xind),bootN);
    sigmall_arr{Xind} = zeros(Ddsub(Xind),bootN);
    psiall_arr{Xind} = zeros(d,bootN);
    phiall_arr{Xind} = zeros(Qdsub(Xind),bootN);
    omegall_arr{Xind} = zeros(Rdsub(Xind),bootN);
    logpall_arr{Xind} = zeros(1,bootN);
    
end

n = sum(Nsub);
Xall = zeros(n,d);

k = 1;
for Xind = 1:Xnm
    Xall(k:k+Nsub(Xind)-1,:) = Xarr{Xind};
    k = k + Nsub(Xind);
end

h = waitbar(0/bootN,'Bootstrapping FA parameters for permutation test');

i = 0;
n_fail = 0;

while i < bootN
    
    bootinds = randperm(n);
    k = 1;
    
    itersub = zeros(Xnm,1);
    logpsub = itersub;
    
    psiarrtmp = {};
    Warrtmp = {};
    sigmarrtmp = {};
    phiarrtmp = {};
    omegarrtmp = {};
    
    for Xind = 1:Xnm
        
        Xtmp = Xall(bootinds(k:k+Nsub(Xind)-1),:);
        k = k + Nsub(Xind);
        [Wtmp,psitmp,phitmp,sigmatmp,omegatmp,itertmp] = ...
            CFAfull_spc(Xtmp,paramstrct_arr{Xind}.Q,100,msub(Xind),paramstrct_arr{Xind}.D,paramstrct_arr{Xind}.R,0);
        
        itersub(Xind) = itertmp;
        
        psiarrtmp{Xind} = psitmp;
        Warrtmp{Xind} = Wtmp;
        sigmarrtmp{Xind} = sigmatmp;
        phiarrtmp{Xind} = phitmp;
        omegarrtmp{Xind} = omegatmp;
        
        if size(paramstrct_arr{Xind}.glob,2) > 1
            Warrtmp{Xind} = rotate_maxaccum(Warrtmp{Xind},paramstrct_arr{Xind}.glob);
        end
        
        Chat = Wtmp*Wtmp' + paramstrct_arr{Xind}.D*diag(sigmatmp)*paramstrct_arr{Xind}.D'...
            + paramstrct_arr{Xind}.Q*diag(phitmp)*paramstrct_arr{Xind}.Q' + ...
            diag(psitmp) + paramstrct_arr{Xind}.R*diag(omegatmp)*paramstrct_arr{Xind}.R';
        
        logpsub(Xind) = -n/2*(d*log(2*pi) + log(det(Chat)) + trace(inv(Chat)*Carr{Xind}));
        
    end
    
    if sum(itersub<100000)==Xnm
        
        i = i + 1;
        
        for Xind = 1:Xnm
            
            if msub(Xind) > 0
                Wall_arr{Xind}(:,:,i) = Warrtmp{Xind};
            end
            
            if Ddsub(Xind) > 0
                sigmall_arr{Xind}(:,i) = sigmarrtmp{Xind};
            end
            
            psiall_arr{Xind}(:,i) = psiarrtmp{Xind};
            
            if Qdsub(Xind) > 0
                phiall_arr{Xind}(:,i) = phiarrtmp{Xind};
            end
            
            if Rdsub(Xind) > 0
                omegall_arr{Xind}(:,i) = omegarrtmp{Xind};
            end
            
            logpall_arr{Xind}(i) = logpsub(Xind);
            
        end
        
        
        
        h = waitbar(i/bootN,h,'Bootstrapping FA parameters for permutation test');
        
    else
        n_fail = n_fail + 1
        
    end
    
end

permstrct.Wall_arr = Wall_arr;
permstrct.sigmall_arr = sigmall_arr;
permstrct.psiall_arr = psiall_arr;
permstrct.phiall_arr = phiall_arr;
permstrct.omegall_arr = omegall_arr;
permstrct.logall_arr = logpall_arr;

close(h)

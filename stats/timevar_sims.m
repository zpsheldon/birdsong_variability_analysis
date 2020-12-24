clear
clc
close all

seqN1 = 50;
seqN2 = 100;
simN = 20;
Mvc = 3:10;

Mvc = 7;

for mind = 1:numel(Mvc)
   
    M = Mvc(mind);
    D = -diff(eye(M+1))';
    Q = zeros(M,0);
    
    Wrlmat = zeros(simN,M);
    psirlmat = Wrlmat;
    Whatmat1 = Wrlmat;
    psihatmat1 = Wrlmat;
    Whatmat2 = Wrlmat;
    psihatmat2 = Wrlmat;
    
    for simind = 1:simN
        
        simind
        lenrl = 30+150*rand(1,M);
        Wrl = (3/200)*lenrl + randn(1,M);
        psirl = diag(5*rand(1,M));
        sigmarl = diag(.5*rand(1,M+1));
        Wrlmat(simind,:) = Wrl;
        psirlmat(simind,:) = diag(psirl);
        
        Xmn = repmat(lenrl,seqN1,1);
        
        z = randn(seqN1,1);
        u = randn(seqN1,M+1);
        eta = randn(seqN1,M);
        
        X = Xmn + z*Wrl + u*sigmarl*D + eta*psirl;
        
        [What,psihat,phihat,sigmahat,omegahat] = CFAfull_prior_spc(X,Q,100,1);
        Whatmat1(simind,:) = What;
        psihatmat1(simind,:) = psihat;
        
        
        Xmn = repmat(lenrl,seqN2,1);
        
        z = randn(seqN2,1);
        u = randn(seqN2,M+1);
        eta = randn(seqN2,M);
        
        X = Xmn + z*Wrl + u*sigmarl*D + eta*psirl;
        
        [What,psihat,phihat,sigmahat,omegahat] = CFAfull_prior_spc(X,Q,100,1);
        Whatmat2(simind,:) = What;
        psihatmat2(simind,:) = psihat;
        %       [What,psihat,phihat,sigmahat,omegahat] = CFAfull_spc(X,Q,100,1);
        

        
%        indstmp = find(abs(sqrt(psihat)-1)<.1 & diag(psirl)>2.5);
%         
%         if ~isempty(indstmp)
%            'here' 
%         end
        
    end
    
    
end

figure; 
subplot(2,2,1)
plot(Wrlmat(:),Whatmat1(:),'o'); hold on; plot([0 4],[0 4],'k--')

subplot(2,2,2)
plot(vec(psirlmat(:,2:end-1)),sqrt(vec(psihatmat1(:,2:end-1))),'o'); 
hold on
plot(vec(psirlmat(:,[1 end])),sqrt(vec(psihatmat1(:,[1 end]))),'o'); 
hold on; plot([0 5],[0 5],'k--')


subplot(2,2,3)
plot(Whatmat1(:),Whatmat2(:),'o'); hold on; plot([0 4],[0 4],'k--')

subplot(2,2,4)
plot(sqrt(vec(psihatmat1(:,2:end-1))),sqrt(vec(psihatmat2(:,2:end-1))),'o'); 
hold on
plot(sqrt(vec(psihatmat1(:,[1 end]))),sqrt(vec(psihatmat2(:,[1 end]))),'o'); 
hold on; plot([0 5],[0 5],'k--')
function [h,Wmat,psimat,sigmamat,lenmat] = plot_prepost4(seqarr,lenarr,xvc)

h = zeros(1,5);

clrarr = {'k','b','r','m','g'};

if nargin == 2
    xvc = 1:length(lenarr);
end

Wmat = zeros(length(seqarr),length(lenarr));
psimat = Wmat;
sigmamat = Wmat;
lenmat = Wmat;

M = length(seqarr);

Q = mkQmat(seqarr);

for i = 1:length(lenarr)
    
    N = size(lenarr{i},1);
    
    lens = lenarr{i};
    
    if size(lens,1) > 20
        
        [W,psi,phi,sigma,iter,n_fail,logp] = CFAJ_spc(lens,Q,100);
        
        maxind = find(logp == max(logp));
        Wmat(:,i) = W(:,maxind);
        psimat(:,i) = sqrt(psi(:,maxind));
        sigmamat(:,i) = sqrt(sigma(:,maxind));
        lenmat(:,i) = mean(lens)';
        
    end
    
end

figure
h(3) = axes;
clrarr = {'k','b','r','m','g','c'};
markarr = {'o','+','x','d','s'};

lgndtxt = seqarr;
lgndtxt{end+1} = 'mean';

hold on
for i = 1:size(lenmat,1)
    clrind = modspec(i,length(clrarr));
    markind = modspec(i,length(markarr));
    plot(xvc,Wmat(i,:),[clrarr{clrind} '-' markarr{markind}],'MarkerFaceColor',clrarr{clrind});
end

% errorbar(seqlen,seqlenerr,'k','linewidth',3)
plot(xvc,mean(Wmat),'k-','linewidth',3)
legend(lgndtxt)


title('tempo noise')

figure
h(4) = axes;
clrarr = {'k','b','r','m','g','c'};
markarr = {'o','+','x','d','s'};

lgndtxt = seqarr;
lgndtxt{end+1} = 'mean';

hold on
for i = 1:size(lenmat,1)
    clrind = modspec(i,length(clrarr));
    markind = modspec(i,length(markarr));
    plot(xvc,psimat(i,:),[clrarr{clrind} '-' markarr{markind}],'MarkerFaceColor',clrarr{clrind});
end

% errorbar(seqlen,seqlenerr,'k','linewidth',3)
plot(xvc,mean(psimat),'k-','linewidth',3)
legend(lgndtxt)

title('fast noise')



figure
h(5) = axes;
clrarr = {'k','b','r','m','g','c'};
markarr = {'o','+','x','d','s'};

lgndtxt = seqarr;
lgndtxt{end+1} = 'mean';

hold on
for i = 1:size(lenmat,1)
    clrind = modspec(i,length(clrarr));
    markind = modspec(i,length(markarr));
    plot(xvc,sigmamat(i,:),[clrarr{clrind} '-' markarr{markind}],'MarkerFaceColor',clrarr{clrind});
end

% errorbar(seqlen,seqlenerr,'k','linewidth',3)
plot(xvc,mean(sigmamat),'k-','linewidth',3)
legend(lgndtxt)


title('jitter')
function plot_prepost(seqarr,lenspre,lenspost,lenspost2)

M = length(seqarr);
ci = zeros(2,M);
for i = 1:M
    [h,p,ci(:,i)] = ttest2(lenspost(:,i),lenspre(:,i));
end

figure
subplot(2,1,1)

d = mean(lenspost)-mean(lenspre);
errorbar(1:M,d,d-ci(1,:),ci(2,:)-d,'b')
hold on

if nargin == 4
    ci2 = zeros(2,M);
    for i = 1:M
        [h,p,ci2(:,i)] = ttest2(lenspost2(:,i),lenspost(:,i));
    end
    d2 = mean(lenspost2)-mean(lenspost);
    errorbar(1:M,d2,d2-ci2(1,:),ci2(2,:)-d2,'r')
    
    legend('post1','post2')
end

plot([1 length(seqarr)],[0 0],'k--')
set(gca,'xtick',1:length(seqarr),'xticklabel',seqarr)
ylabel('post - pre mean len (msec)')


Npre = size(lenspre,1);
Npost = size(lenspost,1);


CL_pre = sqrt((Npre-1)*var(lenspre)/chi2inv(.025,(Npre-1))) ./ mean(lenspre); 
CU_pre = sqrt((Npre-1)*var(lenspre)/chi2inv(.975,(Npre-1))) ./ mean(lenspre); 

CL_post = sqrt((Npost-1)*var(lenspost)/chi2inv(.025,(Npost-1))) ./ mean(lenspost); 
CU_post = sqrt((Npost-1)*var(lenspost)/chi2inv(.975,(Npost-1))) ./ mean(lenspost); 


subplot(2,1,2)
errorbar(1:M,cv(lenspre),cv(lenspre)-CL_pre,CU_pre-cv(lenspre),'k')
hold on
errorbar(1:M,cv(lenspost),cv(lenspost)-CL_post,CU_post-cv(lenspost),'b')

if nargin == 4
    Npost2 = size(lenspost2,1);   
    CL_post2 = sqrt((Npost2-1)*var(lenspost2)/chi2inv(.025,(Npost2-1))) ./ mean(lenspost2);
    CU_post2 = sqrt((Npost2-1)*var(lenspost2)/chi2inv(.975,(Npost2-1))) ./ mean(lenspost2);
    
    errorbar(1:M,cv(lenspost2),cv(lenspost2)-CL_post2,CU_post2-cv(lenspost2),'r')
    legend({'pre','post','post1'})
else
    legend({'pre','post'})
end

set(gca,'xtick',1:length(seqarr),'xticklabel',seqarr)
ylabel('length CV')
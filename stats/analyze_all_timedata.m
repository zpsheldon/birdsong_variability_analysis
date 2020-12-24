clear
clc

psiexp = [];
psinorm = [];
lenexp = [];
lenorm = [];
Wnorm = [];
Wexp = [];
bird = [];
syl = [];
cond = [];

birdspec = [];
sylspec = [];
stdnorm = [];
stdexp = [];
condspec = [];

load('timeanl_dir.mat')
psiexp = [psiexp;paramat(:,2)];
psinorm = [psinorm;paramat(:,1)];
lenorm = [lenorm;paramat(:,3)];
lenexp = [lenexp;paramat(:,4)];
Wnorm = [Wnorm;paramat(:,5)];
Wexp = [Wexp;paramat(:,6)];

syl = [syl;sylall];
bird = [bird;birdall];
cond = [cond;zeros(numel(sylall),1)];

load('specanl_dir.mat')
sylspec = [sylspec;sylall];
birdspec = [birdspec;birdall];
stdexp = [stdexp;paramat(:,2)];
stdnorm = [stdnorm;paramat(:,1)];

condspec = [condspec;zeros(numel(sylall),1)];

load('timeanl_LCstim.mat')


load('timeanl_NE.mat')
psiexp = [psiexp;paramat(:,2)];
psinorm = [psinorm;paramat(:,1)];
lenorm = [lenorm;paramat(:,3)];
lenexp = [lenexp;paramat(:,4)];
syl = [syl;sylall];
bird = [bird;birdall];
cond = [cond;ones(numel(sylall),1)];
Wnorm = [Wnorm;paramat(:,5)];
Wexp = [Wexp;paramat(:,6)];

load('specanl_NE.mat')
sylspec = [sylspec;sylall];
birdspec = [birdspec;birdall];
stdexp = [stdexp;paramat(:,2)];
stdnorm = [stdnorm;paramat(:,1)];
condspec = [condspec;ones(numel(sylall),1)];


load('timeanl_PHE.mat')

d = psiexp-psinorm;

dmat = [];
dmat2 = [];
sylarr = [];
bird2 = [];

dmatspec = [];
sylarrspec = [];
birdspec2 = [];

birdU = unique(bird);

for birdind = 1:numel(birdU)
    
    birdinds0 = strcmp(bird,birdU(birdind)) & cond==0;
    birdinds1 = strcmp(bird,birdU(birdind)) & cond==1;
    
    if sum(birdinds0)>0 & sum(birdinds1)>0
        sylarr = [sylarr;[syl(birdinds0) syl(birdinds1)]];
        dmat = [dmat;psinorm(birdinds0) psiexp(birdinds0) psinorm(birdinds1) psiexp(birdinds1)];
        dmat2 = [dmat2;Wnorm(birdinds0) Wexp(birdinds0) Wnorm(birdinds1) Wexp(birdinds1)];
        bird2 = [bird2;repmat(birdU(birdind),sum(birdinds0),1)];
        
    end
    
    birdinds0 = strcmp(birdspec,birdU(birdind)) & condspec==0;
    birdinds1 = strcmp(birdspec,birdU(birdind)) & condspec==1;
    if sum(birdinds0)>0 & sum(birdinds1)>0
        
        sylarrspec = [sylarrspec;[sylspec(birdinds0) sylspec(birdinds1)]];
        dmatspec = [dmatspec;stdnorm(birdinds0) stdexp(birdinds0) stdnorm(birdinds1) stdexp(birdinds1)];
        birdspec2 = [birdspec2;repmat(birdU(birdind),sum(birdinds0),1)];
        
    end
    
end

d1 = sqrt(dmat(:,2))-sqrt(dmat(:,1));
d2 = sqrt(dmat(:,4))-sqrt(dmat(:,3));
figure;
subplot(2,1,1)
plot(d1,d2,'o')

d1spec = dmatspec(:,2)-dmatspec(:,1);
d2spec = dmatspec(:,4)-dmatspec(:,3);

subplot(2,1,2)
plot(d1spec,d2spec,'o')

dmat2 = [];
% 
% for i = 1:size(dmatspec,2)
%    inds = strcmp(bird,birdspec(i))
%     
% end

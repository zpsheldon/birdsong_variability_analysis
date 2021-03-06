function s = get_seq_timing_amp(sylarr,seqvc,templatestrct,seqstrct,wavdir,N,propstrct,dirstrct)

sylnm = length(sylarr);
oninds = findSeq(seqstrct.vc,seqvc);
seqnm = length(oninds);

tmstmp = seqstrct.cliptms;
lenstmp = seqstrct.cliplens;
wavindstmp = seqstrct.wavinds;

if strcmp(sylarr(1),'*')
    oninds = oninds+1;
    sylnm = sylnm-1;
    sylarr = sylarr(2:end);
end

if strcmp(sylarr(end),'*')
    sylnm = sylnm-1;
    sylarr = sylarr(1:end-1);
end

notenm = 0;
for sylind = 1:length(sylarr)
    if ~strcmp(sylarr{sylind},'sil')
        notenm = notenm + length(templatestrct.featms{find(strcmp(templatestrct.speclabs,sylarr{sylind}))}) - 1;
    end
end

sylvc = zeros(1,notenm);
notenmvc = zeros(1,sylnm);

[b,a] = butter(propstrct.filt_order,2*[propstrct.freqmin propstrct.freqmax]/propstrct.fs);
smthbins = 2^round(log2(propstrct.smthwn / propstrct.timewin_fn));
freqinds = templatestrct.freqinds;
pad = zeros(size(templatestrct.specarr{1},1),smthbins+1);

mfct = 1/2;
minlim = round( propstrct.fs * 15 / (propstrct.f_winadv * 1000));

sampinds = randperm(seqnm);
sampinds = sampinds(1:min(N,seqnm));
sampnm = length(sampinds);


lens = zeros(sampnm,notenm);
ons = zeros(sampnm,sylnm);
sylens = ons;


lastind = 0;

for sylind = 1:sylnm
    
    sylid = sylarr{sylind};
    
    if ~strcmp(sylid,'sil')
        
        template_syl_ind = find(strcmp(templatestrct.speclabs,sylarr{sylind}));
        
        template = templatestrct.specarr{template_syl_ind};
        
        s = zeros(sampnm,size(template,2));
        
        h = waitbar(0/sampnm,['Mapping syl ' sylid]);
        
        for sampind = 1:sampnm
            seqind = sampinds(sampind);
            clipind = oninds(seqind)+sylind-1;
            
            samp1 = ceil(seqstrct.cliptms(clipind) * propstrct.fs / 1000);
            samp2 = samp1 + floor(seqstrct.cliplens(clipind) * propstrct.fs / 1000) - 2;
            flnm = [wavdir filesep dirstrct.wavfls{seqstrct.wavinds(clipind)}];
            
            s_dims = wavread(flnm,'size');
            samp2 = min(samp2,max(s_dims));
            
            signal = filter(b,a,wavread(flnm,double([samp1,samp2])));
            
            spec = tdft(signal,propstrct.f_winlen,propstrct.f_winadv);
            spec = diff(smoothVecGauss([pad spec(freqinds,:) pad],smthbins)')';
            
            mlim = max(minlim,round(size(spec,2)*mfct));
            
            w = dtwProd(sum(template),sum(spec),'pwt',[1.5 1.5 1],'mlim',mlim);
            
            [dmy,winds] = sort(w(:,1));
            w = w(winds,:);
            w = w(w(:,1) > 0,:);
            
            s(sampind,w(:,1)) = w(:,2)';
            
            h = waitbar(sampind/sampnm,h);
            
        end
        
        close(h)
        
        s = s * 1000 * propstrct.f_winadv / propstrct.fs;
        
        lenstmp = chunk(diff(s')',templatestrct.featms{template_syl_ind});
        onstmp =  s(:,templatestrct.featms{template_syl_ind}(1));
        
        noteinds = lastind+1:lastind+length(templatestrct.featms{template_syl_ind})-1;
        lastind = noteinds(end);
        
        lens(:,noteinds) = lenstmp;
        ons(:,sylind) = onstmp + tmstmp(oninds(sampinds)+sylind-1);
        
        if length(noteinds)>1
            sylens(:,sylind) = sum(lens(:,noteinds)')';
        else
            sylens(:,sylind) = lens(:,noteinds);
        end
        
        sylvc(noteinds) = sylind;
        notenmvc(sylind) = length(templatestrct.featms{template_syl_ind})-1;
        
    end
    
end

gaps = diff(ons')';
gaps = gaps - sylens(:,1:end-1);

intnm = size(lens,2)+size(gaps,2);
lenseq = zeros(sampnm,intnm);
sg = ones(1,intnm);
gapinds = cumsum([notenmvc(1:end-1) + ones(1,sylnm-1)]);
sylinds = setdiff(1:intnm,gapinds);
sg(gapinds) = 0;

lenseq(:,sylinds) = lens;
lenseq(:,gapinds) = gaps;

% [lenseq,sngDlt] = cleanInts(lenseq,5,1);

seqarr(1:2:2*sylnm-1) = sylarr;

for i = 2:2:2*sylnm-1
    seqarr{i} = [seqarr{i-1} seqarr{i+1}];
end

notenms = {};
sylind = 1;
noteind = 1;

for i = 1:length(sg)
    if sg(i) == 1
        if notenmvc(find(strcmp(sylarr,seqarr{sylind}))) > 1
            notenms{i} = [seqarr{sylind} '-' num2str(noteind)];
        else
            notenms{i} = seqarr{sylind};
        end
        noteind = noteind + 1;
    else
        sylind = sylind + 1;
        notenms{i} = seqarr{sylind};
        sylind = sylind + 1;
        noteind = 1;
    end
end


clear s

s.sylarr = sylarr;
s.seqarr = notenms;
s.lenseq = lenseq;
s.sg = sg;
s.seqinds = oninds(sampinds);

if exist('sngDlt')
    s.sngDlt = sngDlt;
    s.seqinds = oninds(setdiff(1:length(oninds),sngDlt));
end
function s = get_seq_timing(sylarr,seqvc,featms,seqstrct,sampdir)

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
    notenm = notenm + length(featms{sylind}) - 1; 
end

sylvc = zeros(1,notenm);
notenmvc = zeros(1,sylnm);

lens = zeros(seqnm,notenm);
ons = zeros(seqnm,sylnm);
sylens = ons;

lastind = 0;

for sylind = 1:sylnm
    
    sylid = sylarr{sylind};
    load([sampdir filesep sylid 'warp.mat'],'cliptms','wavinds','cliplens')
    
    flstrct =  whos('-file',[sampdir filesep sylid 'warp.mat']);
    flvars = {flstrct.name};
    if any(strcmp('matnm',flvars))
        load([sampdir filesep sylid 'warp.mat'],'matnm')
    else
        matnm = 1;
    end
        
    
    noteinds = lastind+1:lastind+length(featms{sylind})-1;
    lastind = noteinds(end);
    
    warpinds = zeros(seqnm,1);
    
    for seqind = 1:seqnm
        warpind = find(wavinds == wavindstmp(oninds(seqind)) & cliptms==tmstmp(oninds(seqind)+sylind-1));
        if isempty(warpind)
            errordlg('Warp index not found, bug CG')
            return;
        else
            warpinds(seqind) = warpind;
        end
    end
    
    lenstmp = zeros(length(cliptms),length(featms{sylind})-1);
    onstmp = zeros(length(cliptms),1);
    
    if matnm > 1
        
        ind0 = 1;
        
        for matind = 1:matnm
            sname = ['s' num2str(matind)];
            load([sampdir filesep sylid 'warp.mat'],sname)
            eval(['stmp = ' sname ';'])
            clear(sname)
            
            indstmp = ind0:ind0 + size(stmp,1) - 1;
            ind0 = ind0 + size(stmp,1);
            lenstmp(indstmp,:) = chunk(diff(stmp')',featms{sylind});
            onstmp(indstmp,:) =  stmp(:,featms{sylind}(1));
            
            clear stmp
            
        end
        
    else
         load([sampdir filesep sylid 'warp.mat'],'s')
         lenstmp = chunk(diff(s')',featms{sylind});
         onstmp =  s(:,featms{sylind}(1));
        
         clear s
    end

    lens(:,noteinds) = lenstmp(warpinds,:);
    ons(:,sylind) = onstmp(warpinds) + tmstmp(oninds+sylind-1);
    
    if length(noteinds)>1
        sylens(:,sylind) = sum(lens(:,noteinds)')';
    else
        sylens(:,sylind) = lens(:,noteinds);
    end
    
    sylvc(noteinds) = sylind;
    notenmvc(sylind) = length(featms{sylind})-1;
    
end

gaps = diff(ons')';
gaps = gaps - sylens(:,1:end-1);

intnm = size(lens,2)+size(gaps,2);
lenseq = zeros(seqnm,intnm);
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
s.seqinds = oninds;

if exist('sngDlt')
    s.sngDlt = sngDlt;
    s.seqinds = oninds(setdiff(1:length(oninds),sngDlt));
end
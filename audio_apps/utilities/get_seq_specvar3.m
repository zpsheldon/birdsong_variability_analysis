function [distmat,spmnarr,spstdarr,featms,seqarr,sampinds] = get_seq_specvar3(sylarr,seqvc,featms,seqstrct,wavdir,propstrct,templatestrct,dirstrct,N)

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

[b,a] = butter(propstrct.filt_order,2*[propstrct.freqmin propstrct.freqmax]/propstrct.fs);

sampinds = randperm(seqnm);
sampinds = sampinds(1:min(N,seqnm));
sampnm = length(sampinds);

notenm = 0;
for sylind = 1:length(sylarr)
    notenm = notenm + length(featms{sylind}) - 1;
end

sylvc = zeros(1,notenm);
notenmvc = zeros(1,sylnm);

distmat = zeros(sampnm,notenm);
spmnarr = {};
spstdarr = {};
seqarr = {};

specarrtmp = {};

if sampnm > 0
    
    h = waitbar(0/sampnm,'Loading wav data');
    
    for sampind = 1:sampnm
        
        seqind = sampinds(sampind);
        flnm = [wavdir filesep dirstrct.wavfls{seqstrct.wavinds(oninds(seqind))}];
        
        clipindon = oninds(seqind);
        clipindoff = oninds(seqind)+sylnm-1;
        
        signal_samp1 = ceil(seqstrct.cliptms(clipindon) * propstrct.fs / 1000);
        signal_samp2 = ceil(seqstrct.cliptms(clipindoff) * propstrct.fs / 1000);
        signal_samp2 = signal_samp2 + floor(seqstrct.cliplens(clipindoff) * propstrct.fs / 1000) - 2;
        s_dim = wavread(flnm,'size');
        signal_samp2 = min(signal_samp2,max(s_dim));
        
        signal = wavread(flnm,double([signal_samp1,signal_samp2]));
        signal = filter(b,a,signal);
        
        s_dim = size(signal);
        
        for sylind = 1:sylnm
            
            sylid = sylarr{sylind};
            
            indOn = featms{sylind}(1);
            indOff = featms{sylind}(end);
            
            clipind = oninds(seqind)+sylind-1;
            
            samp1 = ceil(seqstrct.cliptms(clipind) * propstrct.fs / 1000) - signal_samp1 + 1;
            samp2 = samp1 + floor(seqstrct.cliplens(clipind) * propstrct.fs / 1000) - 2;
            samp2 = min(samp2,max(s_dim));
            
            signal_tmp = signal(samp1:samp2);
            
            spec = tdft(signal_tmp,propstrct.f_winlen,propstrct.f_winadv);
            spec = log(max(spec,propstrct.amp_cutoff))-log(propstrct.amp_cutoff);
            specarrtmp{sylind}{sampind} = spec(templatestrct.freqinds,:);
            
            specarrtmp{sylind}{sampind} = specarrtmp{sylind}{sampind} / mean(specarrtmp{sylind}{sampind}(:));
            
        end

        h = waitbar(sampind/sampnm,h,'Loading wav data');
        
    end
    
    close(h)
    
    lastind = 0;
    
    for sylind = 1:sylnm
        [spmn,spstd,spalgn] = mk_spec_template(specarrtmp{sylind},sylarr{sylind});
        
        template = templatestrct.specarr{seqvc(sylind)};
        template = template / mean(template(:));
        
        w=dtwDist(template,spmn,[1.5 1.5 1]);
        
        for featind = 1:length(featms{sylind})
            [dmy,indtmp] = min(abs(w(:,1)-featms{sylind}(featind)));
            featms{sylind}(featind) = round(w(indtmp,2));
        end
        
        indOn = featms{sylind}(1);
        indOff =  featms{sylind}(end);
                
        noteinds = lastind+1:lastind+length(featms{sylind})-1;
        lastind = noteinds(end);
        
        spmn = spmn(:,indOn:indOff);
        spalgn = spalgn(:,:,indOn:indOff);
        
        distmatmp = spalgn - shiftdim(repmat(spmn,[1,1,size(spalgn,1)]),2);
        %     distarr{templateind} = sqrt(sum(sum(distmat.^2,2),3)/numel(spmn));
        %     distarr{templateind} = sum(sum(abs(distmat),2),3)/numel(spmn);
        
        distmp = squeeze((mean(distmatmp.^2,2)));
        
        featms{sylind} = featms{sylind} - featms{sylind}(1) + 1;
        notelentmp = diff(featms{sylind});
        notelentmp(end) = notelentmp(end) + 1;
        
        if length(featms{sylind}) > 2
            distmat(:,noteinds) = chunk(distmp,featms{sylind}) ./ repmat(notelentmp,sampnm,1);
        else
            distmat(:,noteinds) = mean(distmp,2);
        end
        
        spmnarr{sylind} = spmn;
        spstdarr{sylind} = spstd(:,indOn:indOff);
        
        sylvc(noteinds) = sylind;
        notenmvc(sylind) = length(featms{sylind})-1;
        
    end
    
    sylind = 1;
    noteind = 1;
    
    seqarr = {};
    
    for i = 1:notenm
        if notenmvc(sylvc(i)) > 1
            seqarr{i} = [sylarr{sylvc(i)} '-' num2str(noteind)];
        else
            seqarr{i} = sylarr{sylvc(i)};
        end
        
        if noteind==notenmvc(sylvc(i))
            noteind = 1;
        else
            noteind = noteind + 1;
        end
    end
    
end

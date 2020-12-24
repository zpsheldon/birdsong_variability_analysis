function match_directory(wavdir,trainstrct,templatestrct,propstrct,kopt,kthresh)
% hObject    handle to match_all_clips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


flg = check_headers(wavdir);

match_filename = [wavdir 'matchdata.mat'];

load([wavdir filesep 'wavdirinfo.mat']);
% load([wavdir filesep 'clipdirinfo.mat']);

freqinds = trainstrct.freqinds;
labsu = templatestrct.speclabs;
labsu{end+1} = 'x';
% clipnm = length(clipstrct.clipons);

wavnm = length(dirstrct.wavinds);

[b,a] = butter(propstrct.filt_order,2*[propstrct.freqmin propstrct.freqmax]/propstrct.fs);
winlen = propstrct.f_winlen;
winadv = propstrct.f_winadv;

nullen = 60;

dt = 1000 * winadv / propstrct.fs;

nullbins = round(nullen * propstrct.fs / (1000 * winadv));

clipind = 1;

f = propstrct.fs*[0:winlen/2]/winlen;

trainstrct.freqinds = find(f>=propstrct.freqmin & f<=propstrct.freqmax);

f = f(trainstrct.freqinds);

strct = ABconfig;
amp_cutoff = propstrct.amp_cutoff;
speclens = zeros(length(trainstrct.specarr)+1,1);

speclens = trainstrct.speclens;

pcainds = find(~strcmp(trainstrct.speclabs,'x'));

pcaopt = 1;
if pcaopt
    specall = zeros(sum(trainstrct.speclens(pcainds)),size(trainstrct.specarr{1},1));
    k = 1;
    for i = 1:length(pcainds)
        specall(k:k+trainstrct.speclens(pcainds(i))-1,:) = trainstrct.specarr{pcainds(i)}';
        k = k + trainstrct.speclens(pcainds(i));
    end
    
    coeff = princomp(specall);
%     cumvar = cumsum(latent)/sum(latent);
%        
% %     pcadim = max(find(cumvar<=.9)); 

    pcadim = 5;
    
    coeff = coeff(:,1:pcadim);
    
    specmn = mean(specall);
    
    for i = 1:length(trainstrct.specarr)
        trainstrct.specarr{i} = coeff'*(trainstrct.specarr{i}- repmat(specmn',1,trainstrct.speclens(i)));
    end
    
end

sylinds = 1:length(trainstrct.speclabs);

protonm = length(sylinds);

trainstrct.threshvc = kthresh;
trainstrct.kopt = kopt;

silthresh = trainstrct.threshvc(1);

specnull = zeros(length(trainstrct.freqinds),nullbins);

if pcaopt
    specnull = specnull - repmat(specmn',1,nullbins);
    specnull = coeff'*specnull;
end

speclens(end+1:end+trainstrct.kopt) = nullbins;

sylinds = [sylinds length(trainstrct.speclabs)+1:length(trainstrct.speclabs)+trainstrct.kopt];
trainstrct.speclabs(end+1:end+trainstrct.kopt) = {'sil'};
silind_samp = length(trainstrct.speclabs)-trainstrct.kopt+1:length(trainstrct.speclabs);


labsu = unique(trainstrct.speclabs(sylinds));
labnm = length(labsu);

% clipindsval = find(~(strcmp(trainstrct.speclabs,'x')));
clipindsval = 1:length(trainstrct.speclabs);

outind = find(strcmp(labsu,'x'));
silind = find(strcmp(labsu,'sil'));

trainstrct.threshvc = [trainstrct.threshvc*ones(1,length(templatestrct.speclabs)) 0 0];


clipind = 1;
cliptms = zeros(100000,1);
cliplens = cliptms;
wavinds = cliptms;

msg = ['Matching wav files'];
h = waitbar(1/wavnm,msg);

for wavind = 1:wavnm

    s = wavread([wavdir filesep dirstrct.wavfls{wavind}]);
    s = filter(b,a,s);
    X = tdft(s,winlen,winadv);
    clear s
    X = log(max(X(trainstrct.freqinds,:),amp_cutoff)) - log(amp_cutoff);
        
    xlen = size(X,2);
    
    if pcaopt
       X = X - repmat(specmn',1,xlen); 
       X = coeff'*X;    
    end

    Y = zeros(protonm+trainstrct.kopt,xlen)+inf;
    
    for sylind = 1:protonm
        
        spectmp = trainstrct.specarr{sylinds(sylind)};
        y = xcorr21(spectmp,X,1);
        %
        %         y = filter(b2,a2,y);
        %
        %         pkinds = localmaxder(-y-min(-y));
        
        Y(sylind,1:length(y)) = y;
        
    end
    
    y = xcorr21(specnull,X,1);
    Y(end-trainstrct.kopt+1:end,1:length(y)) = repmat(y,trainstrct.kopt,1);
    
    clear X
    Y = Y(clipindsval,:);
    
    [U,I] = sort(Y,1,'ascend');
    
    
    clear Y
    
    I = I(1:trainstrct.kopt,:);
    U = U(1:trainstrct.kopt,:);
    classes = trainstrct.speclabs(clipindsval(sylinds(I)));
    lens = speclens(clipindsval(sylinds(I)));
    
    V = zeros(labnm,xlen);
    L = V;
    S = L;
    xinds = find(U(1,:)~=inf);
    
    for labind = 1:labnm
        classlog = strcmp(classes(:,xinds),labsu{labind});
        V(labind,xinds) = sum(classlog);
        L(labind,xinds) = sum(lens(:,xinds).*classlog) ./ V(labind,xinds);
        S(labind,xinds) = -sum(U(:,xinds).*classlog) ./ V(labind,xinds);
        
    end
    
    [maxval,maxinds] = max(V);
    maxlens = L(sub2ind(size(L),maxinds,1:size(V,2)));
    maxscrs = S(sub2ind(size(L),maxinds,1:size(V,2)));
    
    [sil_i,sil_j] = find(ismember(I(1:ceil(silthresh),xinds),silind_samp));
%     
%     indsval = setdiff(find(maxinds(xinds)~=outind & maxinds(xinds)~=silind),unique(sil_j));
%     indsval2 = find(maxval(xinds(indsval))>=trainstrct.threshvc(maxinds(xinds(indsval))));
    %
    %     indsval = indsval(indsval2);
    
%     indsval = setdiff(find(maxval(xinds)>=trainstrct.threshvc(maxinds(xinds))),sil_j);
    
    indsval = find(maxval(xinds)>=trainstrct.threshvc(maxinds(xinds)));
    
    clear U I V L S
    
    c = zeros(1,xlen);
    v = c;
    l = c;
    s = c-inf;
    c(xinds(indsval)) = maxinds(xinds(indsval));
    l(xinds(indsval)) = maxlens(xinds(indsval));
    s(xinds(indsval)) = maxscrs(xinds(indsval));
    
    %     l(find(c==outind)) = 0;
    
    c(find(c==0)) = silind;
    v(xinds(indsval)) = maxval(xinds(indsval));
    
    %     v(find(c==outind)) = 0;
    
    %     s(find(c==outind)) = -10000+s(find(c==outind));
    %
%     inds = find(c~=silind);
    inds = 1:length(c);
    
%     inds = 1:length(c);
   
    %     [tempons,tempoffs,matchids] = template2clips2(dt*inds',l(inds)'*dt,v(inds)');
    %     inds = inds(matchids);
    
    silinds = find(c==silind);
    
    maxmatlen = 1000;
    matnm = ceil(length(inds)/maxmatlen);
    
    tempons = [];
    tempoffs = [];
    matchids =[];
    
    indbeg = 1;
    matid = 1;
    
    while matid <= matnm && indbeg < length(inds)
        
        if matid<matnm
            indendtarg = inds(min(indbeg + maxmatlen -1,length(inds)));
            [dmy,minindind] = min(abs(indendtarg-silinds));
            silindtmp = silinds(minindind);
            indsindsval = find(inds<silindtmp);
            [dmy,maxindsval] = max(inds(indsindsval)-silindtmp);
            indend = inds(indsindsval(maxindsval));
            indend = find(inds==indend);
            
        else
            indend = length(inds);
        end
        
        if indend - indbeg + 1 > 2000
            indend = indbeg + 2000;
        end
        
%         if useclipsampflg
%            tms = dt*inds(indbeg:indend)';
%            lens = l(inds(indbeg:indend))'*dt;
%            scrs = s(inds(indbeg:indend))';
%         else
%             
%         end
        
               
        [temponstmp,tempoffstmp,matchidstmp] = template2clips2(dt*inds(indbeg:indend)',l(inds(indbeg:indend))'*dt,s(inds(indbeg:indend))');
        
        tempons = [tempons;temponstmp];
        tempoffs = [tempoffs;tempoffstmp];
        matchids = [matchids;matchidstmp+indbeg-1];
        
        indbeg = indend+1;
        
        matid = matid + 1;
        
    end
    
    classinds = c(inds(matchids));
    
    indsval = find(classinds~=outind & classinds~=silind);
    
    
    inds = clipind:clipind+length(indsval)-1;
    cliptms(inds) = tempons(indsval);
    wavinds(inds) = wavind*ones(length(inds),1);
    cliplens(inds) = (tempoffs(indsval)-tempons(indsval));
    class(inds) = labsu(classinds(indsval));
    
    clipind = clipind + length(indsval);
    
    h = waitbar(wavind/wavnm,h);
    
end

close(h)

matchstrct.cliplabs = class;
matchstrct.wavdir = wavdir;
matchstrct.wavinds = wavinds(1:clipind-1)';
matchstrct.cliptms = cliptms(1:clipind-1)';
matchstrct.cliplens = cliplens(1:clipind-1)';
matchstrct.trainingfl = trainstrct.clipfl;

novinds = 1:length(cliptms);
save(match_filename,'matchstrct','novinds')

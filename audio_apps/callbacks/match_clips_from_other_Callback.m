function match_clips_from_other_Callback(hObject, eventdata, handles)
% hObject    handle to load_templates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat');
if ~filename
    return;
else
    load([pathname filename])
end

clipstrct_other = clipstrct;
templatestrct_other = templatestrct;

clipothernm = length(clipstrct_other.speclabs);

clipstrct = get(handles.load_clips,'UserData');
templatestrct = get(handles.make_template,'UserData');
propstrct = get(handles.template_axis,'UserData');

clipinds = find(strcmp(clipstrct.speclabs,'x'));
clipnm = length(clipinds);

clipstrct.ids = {};
clipstrct.scr = zeros(1,clipnm);
matchind = 1;

kans = inputdlg('Choose k');
kopt = str2num(kans{1});

kans = inputdlg('Choose threshold');
kthresh = str2num(kans{1});
protodist = zeros(1,clipothernm);

labsu = templatestrct_other.speclabs;
labsu{end+1} = 'x';

labnm = length(labsu);
clipcnt = zeros(labnm,1);


specall = zeros(sum(clipstrct_other.speclens),size(clipstrct_other.specarr{1},1));
k = 1;
for i = 1:length(clipstrct_other.specarr)
    specall(k:k+clipstrct_other.speclens(i)-1,:) = clipstrct_other.specarr{i}';
    k = k + clipstrct_other.speclens(i);
end

[coeff,score,latent] = princomp(specall);
cumvar = cumsum(latent)/sum(latent);

%     pcadim = max(find(cumvar<=.9));

pcadim = 5;

coeff = coeff(:,1:pcadim);

specmn = mean(specall);

k = 1;
for i = 1:length(clipstrct_other.specarr)
    clipstrct_other.specarr{i} = score(k:k+clipstrct_other.speclens(i)-1,1:pcadim)';
    k = k + clipstrct_other.speclens(i);
end

h = waitbar(1/clipnm,'Matching samples');

for clipind = 1:clipnm
    clip = clipstrct.specarr{clipinds(clipind)};
    
    clip = clip - repmat(specmn',1,size(clip,2));
    clip = coeff'*clip;
    
    for clipotherind = 1:clipothernm
        clipother = clipstrct_other.specarr{clipotherind};
        protodist(clipotherind) = dtwscore(clipother,clip);
    end
    
    [dmy,srtinds] = sort(protodist,'ascend');
    srtinds = srtinds(1:kopt);
    
    for labind = 1:labnm
        clipcnt(labind) = length(find(strcmp(clipstrct_other.speclabs(srtinds),labsu{labind})));
    end
    
    [maxnm,winind] = max(clipcnt);
    
    if length(winind)>1
        distmp = protodist(srtinds);
        winindind = 1;
        winscr = mean(distmp(find(strcmp(labsu{winind(winindind)},specstrct_other.speclabs(srtinds)))));
        
        for indtmp = 2:length(winind)
            scrtmp = mean(distmp(find(strcmp(labsu{winind(winindind)},specstrct_other.speclabs(srtinds)))));
            if scrtmp < winscr
                winindind = indtmp;
                winscr = scrtmp;
            end
        end
        
        winind = winind(winindind);
    end
    
    votes = clipcnt(winind);
    
    if ~strcmp(labsu{winind},'x') && votes >= kthresh
        class{clipind} = labsu{winind};
    else
        class{clipind} = 'x';
    end
    
    h = waitbar(clipind/clipnm,h);
    
end

close(h)

clipstrct.speclabs(clipinds) = class;

set(handles.load_clips,'UserData',clipstrct);
plotspecstrct(handles,handles.clip_axis,handles.clip_slider,clipstrct,1000,5,'label');

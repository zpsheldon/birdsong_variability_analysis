function s = get_seqs(sylarr,matchstrct,gapmax)

sylnm = length(sylarr);
sampdir = matchstrct.wavdir;

sylid = sylarr{1};
clipinds = find(strcmp(matchstrct.cliplabs,sylarr{1}));
ons2 = matchstrct.cliptms(clipinds)';
lens2 = matchstrct.cliplens(clipinds)';

seqids = matchstrct.wavinds(clipinds);
seqnm = length(seqids);

seqinds = zeros(seqnm,sylnm);
wavinds = seqinds;
lens = seqinds;
ons = seqinds;
scrs = seqinds;

seqinds(:,1) = 1:length(clipinds);
wavinds(:,1) = seqids;
ons(:,1) = ons2;
lens(:,1) = lens2;
scrs(:,1) = matchstrct.clipscrs(clipinds);


tot = seqnm * (sylnm - 1);

h = waitbar(0/tot,'Loading sequences');
k = 0;

for sylind = 2:sylnm

    sylid = sylarr{sylind};
    clipinds = find(strcmp(matchstrct.cliplabs,sylid));
   
    ons2 = matchstrct.cliptms(clipinds)';
    lens2 = matchstrct.cliplens(clipinds)';

    offsets = ons(:,sylind-1) + lens(:,sylind-1);

    ia = find(ismember(matchstrct.wavinds(clipinds),seqids));
    clipinds = clipinds(ia);
    ons2 = ons2(ia);
    lens2 = lens2(ia);

    for seqind = 1:seqnm

        seqindtmp = find(matchstrct.wavinds(clipinds) == seqids(seqind));
        seqindtmp = seqindtmp(find(matchstrct.cliptms(clipinds(seqindtmp)) > offsets(seqind) - 20 & ...
            matchstrct.cliptms(clipinds(seqindtmp)) - offsets(seqind) < gapmax(sylind-1)));

        if length(seqindtmp) > 1
            [dmy,seqindtmpind] = min(abs(matchstrct.cliptms(clipinds(seqindtmp)) - offsets(seqind)));
            seqindtmp = seqindtmp(seqindtmpind);
        end

        if ~isempty(seqindtmp)
            seqinds(seqind,sylind) = ia(seqindtmp);
            wavinds(seqind,sylind) = matchstrct.wavinds(clipinds(seqindtmp));
            ons(seqind,sylind) = ons2(seqindtmp);
            lens(seqind,sylind) = lens2(seqindtmp);
            scrs(seqind,sylind) = matchstrct.clipscrs(clipinds(seqindtmp));
        end

        k = k + 1;
        h = waitbar(k/tot,h,'Loading sequences');

    end

end

close(h)

gaps = diff(ons')';
gaps = gaps - lens(:,1:end-1);

intnm = size(lens,2)+size(gaps,2);
seqnm = size(lens,1);
lenseq = zeros(seqnm,intnm);

lenseq(:,1:2:end) = lens;
lenseq(:,2:2:end) = gaps;

[i,j] = find(wavinds == 0);
indsval = setdiff(1:seqnm,unique(i));

seqinds = seqinds(indsval,:);
wavinds = wavinds(indsval,:);
lens = lens(indsval,:);
ons = ons(indsval,:);
scrs = scrs(indsval,:);
seqnm = length(indsval);
lenseq = lenseq(indsval,:);


seqarr(1:2:intnm) = sylarr;

for i = 2:2:intnm
    seqarr{i} = [seqarr{i-1} seqarr{i+1}];
end

clear s

s.sylarr = sylarr;
s.seqarr = seqarr;
s.gapmax = gapmax;
s.seqinds = seqinds;
s.wavinds = wavinds;
s.lens = lens;
s.ons = ons;
s.lenseq = lenseq;
s.scrs = scrs;

function [tempons,tempoffs,matchids,matchscrs2] = template2clips2(matchons,matchlens,matchscrs,overlap_tm)
%template2clips(clipons,clipoffs,matchons,matchlens,matchscrs)
%algorithm to match clips to templates using parallel matrices rather than
%serial loops

if nargin<4
    overlap_tm = 70;
end

[rows,cols] = size(matchons);

matchid = repmat([1:rows]',1,cols);
matchid = matchid(:);
matchons = matchons(:);
matchlens = matchlens(:);
matchscrs = matchscrs(:);

% valind = setdiff(1:length(matchons),find(matchons==-1000));
% matchid = matchid(valind);
% matchons = matchons(valind);
% matchlens = matchlens(valind);
% matchscrs = matchscrs(valind);

matchnm = length(matchons);

matchoffs = matchons + matchlens;

%by default, determine all match onsets and offsets according to template
%timing

%two template matches conflict if they overlap by more than [overlap_tm]. in the case of a conflict,
%take the one with the highest score

%4 kinds of overlap: straddle onset, straddle offset, encompassing, embedded 


scrmat = repmat(matchscrs,1,matchnm);
dltscr = scrmat > scrmat';
clear scrmat

onmat = repmat(matchons,1,matchnm);
offmat = repmat(matchoffs,1,matchnm);

dlton = onmat >= onmat';
dltoff = offmat <= offmat';

dlt1 = (onmat <= offmat' - overlap_tm) & dltscr & dlton;
dlt2 = (offmat >= onmat' + overlap_tm) & dltoff & dltscr;
dlt3 = (onmat <= onmat' & offmat >= offmat') & dltscr;
dlt4 = dltscr & dlton & dltoff;

[x,dlt1] = find(dlt1);
[x,dlt2] = find(dlt2);
[x,dlt3] = find(dlt3);
[x,dlt4] = find(dlt4);


% [x,dlt1] = find(onmat >= onmat' & onmat <= offmat' - overlap_tm & scrmat > scrmat');
% 
% [x,dlt2] = find(offmat >= onmat' + overlap_tm & offmat <= offmat' & scrmat > scrmat');
% 
% [x,dlt3] = find(onmat <= onmat' & offmat >= offmat' & scrmat > scrmat');
% 
% [x,dlt4] = find(onmat >= onmat' & offmat <= offmat' & scrmat > scrmat');


dltind = unique([dlt1;dlt2;dlt3;dlt4]);

clear onmat offmat scrmat

indmat = repmat([1:matchnm]',1,matchnm);

indlt = indmat(dltind);
indval = setdiff(1:matchnm,indlt);

tempons = matchons(indval);
tempoffs = matchoffs(indval);
matchids = matchid(indval);
matchscrs2 = matchscrs(indval);

%if a match onset/offset is within [clipdist] msec of clip onset/offset, recaculate
%as such

% clipdist = 15;
% 
% for i = 1:length(tempons)
%    
%     [d,ind] = min(abs(tempons(i) - clipons));
%     if d < clipdist
%         tempons(i) = clipons(ind);
%     end
%     
%     [d,ind] = min(abs(tempoffs(i) - clipoffs));
%     if d < clipdist
%         tempoffs(i) = clipoffs(ind);
%     end
%         
% end

[tempons,i] = sort(tempons);
tempoffs = tempoffs(i);
matchids = matchids(i);
matchscrs2 = matchscrs2(i);
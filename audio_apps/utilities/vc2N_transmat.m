function [transarr,Xarr,Narr,matchnms,Xarru,Narru,Iarr,Parr,Earr,onarr] = vc2N_transmat(vc,N_min,M_min,backopt)

if nargin == 3
    backopt = 0;
end

[transmat,P,matchnms] = vc2transmat(vc);

if backopt
    P = transmat ./ (ones(size(transmat,1),1)*sum(transmat));
end

[x1,x2] = find(transmat > N_min);
transnm = length(x1);
X = [x1 x2];

Xarr{1} = X;
Parr{1} = P(sub2ind(size(P),x1,x2));
E = -sum(P' .* log2(max(P',.0001)))';
Earr{1} = E(x1);

Xarru = {};
Narru = {};

Nvc = diag(transmat(x1,x2));
Narr{1} = Nvc;
transarr{1} = transmat;
Iarr{1} = x1;
vclen = length(vc);

M = 2;

while transnm >= 1
    
    transmat2 = transmat;
    transmat = zeros(transnm,size(transmat,2));
    onarrtmp = {};
    
     for chnkind = 1:transnm
        oninds = findSeq(vc,X(chnkind,:));

        if backopt
            oninds = oninds(oninds >1);
            ids = vc(oninds-1);
        else
            oninds = oninds(oninds < vclen-M+1);
            ids = vc(oninds+M);     
        end
        
        
        t = tabulate(ids);
        itmp2 = find(t(:,2));
        
        transmat(chnkind,t(itmp2,1)) = transmat(chnkind,t(itmp2,1)) + t(itmp2,2)';
        onarrtmp{chnkind} = oninds;
    end
    
%  
%     
%     if backopt
%         P = transmat ./ (ones(size(transmat,1),1)*sum(transmat));
%     else        
%         P = transmat ./ (sum(transmat')'*ones(1,size(transmat,2)));
%     end
    
    P = transmat ./ (sum(transmat')'*ones(1,size(transmat,2)));
    
    E = -sum(P' .* log2(max(P',.0001)))';
    
    [x1,x2] = find(transmat >= N_min);
    transnm = length(x1);
    
    indsprev = setdiff(1:size(Xarr{M-1},1),x1);
    Xarru{M-1} = Xarr{M-1}(indsprev,:);
    Narru{M-1} = Narr{M-1}(indsprev,:);
    
    X2 = zeros(transnm,M+1);
    
    if backopt
        for i = 1:transnm
            X2(i,:) = [x2(i) X(x1(i),1:end)];
        end     
                
    else
        for i = 1:transnm
            X2(i,:) = [X(x1(i),1:end) x2(i)];
        end
    end
    
    X = X2;
    
    Nvc = diag(transmat(x1,x2));
    
    Xarr{M} = X2;
    Narr{M} = Nvc;
    Iarr{M} = x1;
    
    Parr{M} = P(sub2ind(size(P),x1,x2));
    Earr{M} = E(x1);
    transarr{M} = transmat;
    onarr{M-1} = onarrtmp;
    
    M = M + 1;
    
end

Xarr = Xarr(1:end-1);
Narr = Narr(1:end-1);
Iarr = Iarr(1:end-1);
Parr = Parr(1:end-1);
Earr = Earr(1:end-1);
transarr = transarr(1:end-1);

if nargin > 2 & M >= M_min
    Xarru = Xarru(M_min:end);
    Narru = Narru(M_min:end);
end
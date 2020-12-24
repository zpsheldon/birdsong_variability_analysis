function [Xarr,Narr,Iarr,Parr,Earr,onarr,transarr] = vc2N_transmat_2seq(vc,seq,Mmax,Nmin,backopt)

if nargin == 4
    backopt = 0;
end

onarr = {};

[transmat,P,matchnms] = vc2transmat(vc);

transarr{1} = transmat;

if backopt
    P = transmat ./ (ones(size(transmat,1),1)*sum(transmat));
end

if length(seq) > 1
    if backopt
        x1 = find(matchnms==seq(end-1));
        x2 = find(matchnms==seq(end));
    else
        x1 = find(matchnms==seq(1));
        x2 = find(matchnms==seq(2));
    end

else

    xtmp = find(matchnms==seq(1));

    [x1,x2] = find(transmat >= Nmin);

    if backopt
        indsval = find(x2==xtmp);
    else
        indsval = find(x1==xtmp);
    end
    x1 = x1(indsval);
    x2 = x2(indsval);

end

X = [x1 x2];

transnm = length(x1);

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

indsval_old = 1:transnm;

vc = vc';
vclen = size(vc,2);
seqlen = length(seq);

M = 2;

while transnm >= 1 && M <= Mmax + length(seq)

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

        if ~isempty(oninds)

            t = tabulate(ids);
            itmp2 = find(t(:,2));

            transmat(chnkind,t(itmp2,1)) = transmat(chnkind,t(itmp2,1)) + t(itmp2,2)';
            onarrtmp{chnkind} = oninds;

        end
    end

    transarr{M} = transmat;

    if backopt
        P = transmat ./ max((ones(size(transmat,1),1)*sum(transmat)),.01);
    else
        P = transmat ./ max((sum(transmat')'*ones(1,size(transmat,2))),.01);
    end

    E = -sum(P' .* log2(max(P',.0001)))';

    [x1,x2] = find(transmat >= Nmin);

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

    if M < length(seq)

        if backopt
            indsval = find(sum(abs(X2-repmat(seq(end-M:end),size(X2,1),1))')==0);
        else
            indsval = find(sum(abs(X2-repmat(seq(1:M+1),size(X2,1),1))')==0);
        end

    elseif length(seq) > 1

        if backopt
            indsval = find(sum(abs(X2(:,end-length(seq)+1:end)-repmat(seq,size(X2,1),1))')==0);
        else
            indsval = find(sum(abs(X2(:,1:length(seq))-repmat(seq,size(X2,1),1))')==0);
        end
        
    else
        
        indsval = 1:size(X2,1);

    end

    Xarr{M} = X2(indsval,:);
    Narr{M} = Nvc(indsval);
    Iarr{M} = x1(indsval);
    ptmp = P(sub2ind(size(P),x1(indsval),x2(indsval)));
    etmp = E(x1(indsval));
    Parr{M} = ptmp(:);
    Earr{M} = etmp(:);

    if ~isempty(onarrtmp)
        onarr{M-1} = onarrtmp(indsval_old);
    end

    indsval_old = indsval;


    M = M + 1;

end

Xarr = Xarr(1:end-1);
Narr = Narr(1:end-1);
Iarr = Iarr(1:end-1);
Parr = Parr(1:end-1);
Earr = Earr(1:end-1);
transarr = transarr(1:end-1);
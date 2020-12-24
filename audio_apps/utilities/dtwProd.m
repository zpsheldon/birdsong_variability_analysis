function [w,D,d]=dtwProd(x,y,varargin)
% function [w,D,d]=dtwProd(x,y,varargin)
% DTW function using outerproduct

[rows,N]=size(x);
[rows,M]=size(y);

fct = 1/4;
mlim_def = round(M*fct);


vars = {'pwt','metric','mlim'};
% vardefs = {[1.5 1.5 1],'prod',mlim_def};
% vardefs = {[2 2 1],'prod',mlim_def};
% vardefs = {[1.4 1.4 1],'prod',mlim_def};
vardefs = {[1 1.5 1.5],'prod',mlim_def};

for i = 1:length(vars)
    
    varind = find(strcmp(varargin,vars{i}));
    if ~isempty(varind)
        varval = varargin{varind+1};
    else
        varval = vardefs{i};
    end
    
    eval([vars{i} ' = varval;'])
    
end

switch metric
    case 'prod'
        d = x'*y;
    case 'specorr'
        d = specorr(x,y);
end

%paramters restrict path to inner-third map diagonals (Sakoe-Chiba band)

%cumulative score matrix
D = mkCumProd(d,pwt,M,N,mlim);

[i,j] = find(D);
inds = find(i==N | j ==M);
Dvals = D(sub2ind(size(D),i,j));
[dmy,maxind] = max(Dvals(inds));
n=i(inds(maxind));
m=j(inds(maxind));


k=1;

%backwards path tracing, stored in w

w = zeros(N+M,2);
w(1,:)=[n,m];

while n > 2 & m > 2
    
    P1 = .5*d(n,m) + .25*d(n-1,m) + .25*d(n-1,m-1);
    P2 = .5*d(n,m) + .25*d(n,m-1) + .25*d(n-1,m-1);
    P3 = d(n,m);

    pts = [D(n-2,m-1),D(n-1,m-2),D(n-1,m-1)];
    pthset = pts + [P1,P2,P3] .* pwt;

    [values,number] = max(pthset);
    
    k=k+1;

    switch number
        case 1
            n=n-1;
            m=m-.5;
            w(k,:) = [n,m];
            
            k=k+1;
            n=n-1;
            m=m-.5;
            w(k,:) = [n,m];

        case 2
            n=n-1;
            m=m-2;
            w(k,:) = [n,m];

        case 3
            n=n-1;
            m=m-1;
            w(k,:) = [n,m];

    end

end


% calculate endings. max/min functions not used above in order to save
% processing time.

if w(k,1)>1 & w(k,2)>1

    if w(k,1)==2 & w(k,2)==2
        k = k+1;
        w(k,:) = [n-1,m-1];

    elseif w(k,1)==2

        P2 = .5*d(n,m) + .25*d(n,m-1) + .25*d(n-1,m-1);
        P3 = d(n,m);

        pts = [D(n-1,m-2),D(n-1,m-1)];
        pthset = pts + [P2,P3] .* pwt(2:3);

        [values,number] = max(pthset);

        k=k+1;

        switch number

            case 1
                n=n-1;
                m=m-2;
                w(k,:) = [n,m];

            case 2
                n=n-1;
                m=m-1;
                w(k,:) = [n,m];
        end

    else

        P1 = .5*d(n,m) + .25*d(n-1,m) + .25*d(n-1,m-1);
        P3 = d(n,m);

        pts = [D(n-2,m-1),D(n-1,m-1)];
        pthset = pts + [P1,P3] .* pwt([1 3]);

        [values,number] = max(pthset);

        k=k+1;

        switch number
            case 1
                n=n-1;
                m=m-.5;
                w(k,:) = [n,m];

                k=k+1;
                n=n-1;
                m=m-.5;
                w(k,:) = [n,m];

            case 2
                n=n-1;
                m=m-1;
                w(k,:) = [n,m];

        end

    end

end

w = w(1:k,:);
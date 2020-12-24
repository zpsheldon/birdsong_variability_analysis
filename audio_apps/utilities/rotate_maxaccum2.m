function V = rotate_maxaccum2(W,Vprev)

[D,M] = size(W);
one = ones(D,1);

T1 = W' * one / norm(W' * one);

T = [T1 null(T1')];
V = W*T;

% V1 = W * W' * one / norm(W' * one);
% 
% [U,E] = eig(W*W' - V1*V1');
% [e,ord] = sort(diag(E),'descend');
% U = U(:,ord(1:M-1));
% E = diag(e(1:M-1));
% 
% V = [V1 U*sqrt(E)];

if nargin==2
    [dmy,V2] = procrustes(Vprev(:,2:end),V(:,2:end));
    
    [L, D, M] = svd(Vprev(:,2:end)' * V(:,2:end));
    T = M * L';
    V2 = V(:,2:end)*T;
    
    V = [V1 V2];
end

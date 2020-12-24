function y = vech(X)

i = find(tril(ones(length(X))));
y = X(i);
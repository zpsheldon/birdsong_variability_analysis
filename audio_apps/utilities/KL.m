function D = KL(mu_f,sigma_f,mu_g,sigma_g)

mu_f = mu_f(:);
mu_g = mu_g(:);

d = length(mu_f);

D = .5*(log(det(sigma_f)/det(sigma_g)) + trace(inv(sigma_g)*sigma_f) - d + ...
    (mu_f - mu_g)' * inv(sigma_g) * (mu_f - mu_g));

function D = KL_gaussmix(pi_f,mu_f,sigma_f,pi_g,mu_g,sigma_g)

[k_f,d] = size(mu_f);
k_g = size(mu_g,1);

D = 0;

for kind_f = 1:k_f
    
    mu_f_a = mu_f(kind_f,:);
    sigma_f_a = squeeze(sigma_f(:,:,kind_f));
    
    expD_ff = 0;
    expD_fg = 0;
    
    for kind_f2 = 1:k_f
        mu_f_a2 = mu_f(kind_f2,:);
        sigma_f_a2 = squeeze(sigma_f(:,:,kind_f2));
        expD_ff = expD_ff + pi_f(kind_f2)*exp(-KL(mu_f_a,sigma_f_a,mu_f_a2,sigma_f_a2));
    end
    
    for kind_g = 1:k_g
        mu_g_b = mu_g(kind_g,:);
        sigma_g_b = squeeze(sigma_g(:,:,kind_g));
        expD_fg = expD_fg + pi_g(kind_g)*exp(-KL(mu_f_a,sigma_f_a,mu_g_b,sigma_g_b));
    end
    
    D = D + pi_f(kind_f) * log(expD_ff / expD_fg);
    
end
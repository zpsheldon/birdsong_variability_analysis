function [X,z,u,eta] = timing_mcsim(N,W,sigma,psi)
% X = timing_mcsim(N,W,sigma,psi)
%
% MATLAB S4: Runs Monte Carlo simulation to geneterate matrix of interval lengths according to
% Glaze & Troyer (2012) timing variability model. 
%
% N is number of sequences generated and W, sigma and psi are parameters of
% the model (see other functions or manuscript for variable definitions).
%
% Output X is N by length(psi) matrix of deviations from the mean length for length(psi) intervals
%
%
% AUTHORS: Chris Glaze (cmglaze@gmail.com) and Todd Troyer (todd.troyer@utsa.edu)

% make sigma and psi into vectors if not already
if min(size(sigma))>1
    sigma = diag(sigma);
end
if min(size(psi))>1
	psi = diag(psi);
end

K = length(psi);
D = -diff(eye(K))';
[Ktmp,Q] = size(W);
Ktmp2 = length(sigma)+1;

if any(abs(diff([K Ktmp Ktmp2])) > 0)
    error('Parameter dimensions do not match. Aborting.'); 
end

% generate latent variables and data
z = randn(N,Q);
u = randn(N,K-1)*diag(sqrt(sigma));
eta = randn(N,K)*diag(sqrt(psi));
X = z*W' + u*D' + eta;

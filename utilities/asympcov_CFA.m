function [C_theta,semW,semN,semJ,semWvc,semNvc,semJvc] = asympcov_CFA(W,psi,sigma,N)

d = length(W);

D = -diff(eye(d))';

theta_sigma = sqrt(diag(sigma));
theta_psi = sqrt(diag(psi));

% dtheta_sigma = (kron(eye(d-1),2*theta_sigma') .* (ones(d-1,1)*vec(eye(d-1))'))*kron(D',D');
% dtheta_psi = kron(eye(d),2*theta_psi') .* (ones(d,1)*vec(eye(d))');

dtheta_sigma = (kron(eye(d-1),2*theta_sigma') .* (ones(d-1,1)*vec(eye(d-1))'))*kron(D',D');
dtheta_psi = kron(eye(d),2*theta_psi') .* (ones(d,1)*vec(eye(d))');


dtheta_W = kron(eye(d),W') + kron(W',eye(d));

dtheta = [dtheta_W;dtheta_sigma;dtheta_psi];


Chat = W*W' + D*sigma*D' + psi;
C_theta = 2*inv(dtheta*inv(kron(Chat,Chat))*dtheta')/N;

sems = sqrt(diag(C_theta));

semWvc = sems(1:d);
semJvc = sems(d+1:2*d-1);
semNvc = sems(2*d:end);

semW = mean(sems(2:d-1));
semJ = mean(sems(d+2:2*d-2));
semN = mean(sems(2*d+1:end-1));

% I = eye(d);
% 
% C_W = 2*(Einv*(W'*Einv*W) + kron(Einv*W,W'*Einv));
% 
% C_W_D = dE_W * kron(Einv,Einv) * [dE_sigma]';
% C_W_I = dE_W * kron(Einv,Einv) * [dE_psi]';

function [Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta)

%Construct coefficient matrices
bbeta = theta(1);
ggamma = theta(2);
llambda = theta(3);
piStar = theta(4);
zetaP = theta(5);
nu = theta(6);
rhoPhi = theta(7);
rhoLambda = theta(8);
rhoZ = theta(9);
sigmaPhi = theta(10);
sigmaLambda = theta(11);
sigmaZ = theta(12);
sigmaR = theta(13);


%%%%%%%%%%%%%% Helper variables
kappa_P = (1-zetaP*bbeta)*(1-zetaP)/zetaP;
psi_P = (1+kappa_P/(bbeta*(1+nu)))^(-1);
lsh = 1/(1+llambda);


%%%%%%%%%%%% Model Solution in terms of coefficient matrices
Phi1 = diag([rhoPhi, rhoLambda, rhoZ, 0,0]);

Phi1(5,:) = [-(kappa_P*psi_P/bbeta)/(1-psi_P*rhoPhi); ...
    - (kappa_P*psi_P/bbeta)/(1-psi_P*rhoLambda);...
    (rhoZ*psi_P)/(1-psi_P*rhoZ); ...
    - psi_P*sigmaR; 0]';

Phi_eps = diag([ sigmaPhi, sigmaLambda, sigmaZ, 1]);
Phi_eps(5,:) = 0;


Psi0 = [log(ggamma),log(lsh),log(piStar),log(piStar*ggamma/bbeta)];

Psi1_X = Phi1(5,:) + [0,0,1,0,-1];

Psi1_lsh = [(1 - ( (1+nu)*kappa_P*psi_P/bbeta)/(1-psi_P*rhoPhi)),...
    - ( (1+nu)*kappa_P*psi_P/bbeta)/(1-psi_P*rhoLambda),...
    + (1+nu)*rhoZ*psi_P/(1-psi_P*rhoZ),...
    - (1+nu)*psi_P*sigmaR, 0];
Psi1_pi = [(1 - ( (1+nu)*kappa_P*psi_P/bbeta)/(1-psi_P*rhoPhi))*(kappa_P/(1-bbeta*rhoPhi)),...
    (1 - ( (1+nu)*kappa_P*psi_P/bbeta)/(1-psi_P*rhoLambda))*(kappa_P/(1-bbeta*rhoLambda)),...
    kappa_P*(1+nu)*rhoZ*psi_P/( (1-psi_P*rhoZ)*(1-bbeta*rhoZ)),...
    -kappa_P*(1+nu)*psi_P*sigmaR, 0];
Psi1_R = [1/bbeta*( (1 - ( (1+nu)*kappa_P*psi_P/bbeta)/(1-psi_P*rhoPhi))*(kappa_P/(1-bbeta*rhoPhi))),...
    1/bbeta*(1 - ( (1+nu)*kappa_P*psi_P/bbeta)/(1-psi_P*rhoLambda))*(kappa_P/(1-bbeta*rhoLambda)),...
    kappa_P*(1+nu)*rhoZ*psi_P/(bbeta*(1-psi_P*rhoZ)*(1-bbeta*rhoZ)),...
    (1-kappa_P*(1+nu)*psi_P/bbeta)*sigmaR,0];


Psi1 = [Psi1_X; Psi1_lsh; Psi1_pi; Psi1_R];

end

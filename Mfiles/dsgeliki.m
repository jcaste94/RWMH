function [liki] = dsgeliki(para)

% This function computes the likelihood of the DSGE model.
% Input = para: vector of structural parameters
% Output= likelihood function

load Data/Y_sim

%reverse transformations in para:
para(1) = 1/(1 + para(1)/100 );
para(2) = exp(para(2)/100);
para(4) = exp(para(4)/100);
para(6) = 1/para(6) - 1;
para(10:13) = para(10:13)/100;

% Solve DSGE model
[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(para);

% no measurement errors
Sigma_u = zeros(4,4);

% evaluate likelihood function
liki = kalman(Psi0 ,Psi1, zeros(4), eye(5,5) , Phi_eps*Phi_eps',Phi1,Y_sim);
liki = sum(liki);

end




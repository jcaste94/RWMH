function [Y, s_t] = DSGE_simulate(theta, Tlength, Tburnin)
% code to simulate observations from a DSGE model
% input is a vector of parameters, and output is a matrix of simulated
% observables

% ----------------
% 1. Housekeeping
% ----------------

[Phi1, Phi_eps, Psi0, Psi1] =  DSGE_soln_matrices(theta);

rhoPhi = theta(7);
rhoLambda = theta(8);
rhoZ = theta(9);
sigmaPhi = theta(10);
sigmaLambda = theta(11);
sigmaZ = theta(12);

% ----------------
% 2. Simulate s_t
% ----------------

% 2.1. Initial conditions
P0 = diag([sigmaPhi^2/(1-rhoPhi^2), sigmaLambda^2/(1-rhoLambda^2), sigmaZ^2/(1-rhoZ^2), 1]);

Phi_past = Phi1(5,1:4); %extract the coefficients for x_{t-1}

s0 = mvnrnd(zeros(4,1),P0)';
x0 = Phi_past*s0;

s0 = mvnrnd(zeros(4,1),P0)'; % now draw the four states again to use as the initial values for the states

s0 = [s0; x0];


% 2.2. Transition equation
eps_t = randn(Tburnin+Tlength,4);
s_t = zeros(Tburnin+Tlength,5);
s_t(1,:) = s0;
for t = 2:(Tlength+Tburnin)
    s_t(t,:) = Phi1*s_t(t-1,:)' + Phi_eps*eps_t(t-1,:)';
end

s_t = s_t(length(s_t)-Tlength+1:length(s_t),:); % throw away burning sample

% ---------------------------
% 3. Turn S into observables
% ---------------------------

% 3.1. Measurement equation
Y = zeros(Tlength,4);
for t = 1:Tlength
    Y(t,:) = Psi0 + ( Psi1*s_t(t,:)' )';
end

end
%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%        Constructing the Candidate Density for MH Algorithm
%
%
%
% Author: Luigi Bocola         lbocola@sas.upenn.edu
% Date  : 06/16/2013
%
% Edited: Juan Castellanos Silván
% Date  : 02/27/2020
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================

clear all
clc
close all
delete *.asv

tic

l = path;
path('Mfiles',path);
path('Optimization Routines',path);
path('Matfiles',path);

disp('                                                                  ');
disp('    BAYESIAN ESTIMATION OF DSGE MODEL: THE CANDIDATE DENSITY      ');
disp('                                                                  ');

%=========================================================================
%                  True parameter values for simulation
%=========================================================================

bbeta    = 1/1.01;
ggamma   = exp(.005);
llambda  = .15;
piStar   = exp(.005);
zetaP    = .65;
nu       = 0;
rhoPhi   = .94;
rhoLambda= .88;
rhoZ     = .13;
sigmaPhi = .01;
sigmaLambda = .01;
sigmaZ   = .01;
sigmaR   = .01;

%define theta vector
theta_true(1)  = bbeta;
theta_true(2)  = ggamma;
theta_true(3)  = llambda;
theta_true(4)  = piStar;
theta_true(5)  = zetaP;
theta_true(6)  = nu;
theta_true(7)  = rhoPhi;
theta_true(8)  = rhoLambda;
theta_true(9)  = rhoZ;
theta_true(10) = sigmaPhi;
theta_true(11) = sigmaLambda;
theta_true(12) = sigmaZ;
theta_true(13) = sigmaR;

% simulate data
[Y_sim, s_sim] = DSGE_simulate(theta_true, 500, 100);

save('Data/Y_sim','Y_sim')

% transform parameters for estimation:
%  1.) 100(1/Beta-1) 
%  2.) 100*log(Gamma)
%  3.) Lambda 
%  4.) 100*log(PiStar) 
%  5.) zetaP 
%  6.) 1/(1+nu) 
%  7.) rhoPhi 
%  8.) rhoLambda 
%  9.) rhoZ 
%  10.) 100*sigmaPhi 
%  11.) 100*sigmaLambda 
%  12.) 100*sigmaZ 
%  13.) 100*sigmaR 

param = theta_true;

param(1)     = 100*(1/param(1) - 1);
param(2)     = 100*log(param(2));
param(4)     = 100*log(param(4));
param(6)     = 1/(1+param(6));
param(10:13) = 100*param(10:13);

disp('                                                                  ');
disp('         *******STEP 1: RECOVERING THE POSTERIOR MODE....*********');
disp('                                                                  ');
 
% we are minimizing the objective function, which is why we are multiplying
% by (-1)
objective = @(theta) prior(theta)*(-1) + dsgeliki(theta)*(-1);
[fh,param_mode,gh,H,itct,fcount,retcodeh] = csminwel(objective,param,eye(length(param)),[] ,10^(-5),200);

% extract the argmin
theta_mode = param_mode;

% un-do the transformations
theta_mode(1)     = 1/(1 + theta_mode(1)/100 );
theta_mode(2)     = exp(theta_mode(2)/100);
theta_mode(4)     = exp(theta_mode(4)/100);
theta_mode(6)     = 1/theta_mode(6) - 1;
theta_mode(10:13) = theta_mode(10:13)/100;

disp('                                                                  ');    
disp('                            THE TRUE VALUES ARE:                  ');
disp('                                                                  ');
disp('   BETA      GAMMA      LAMBDA        Pi*        ZETA_p          NU    RHO_phi      RHO_lambda    RHO_Z  SIGMA_Phi  SIGMA_lambda  SIGMA_Z     SIGMA_R ');
disp(num2str(theta_true))
disp('                                                                  ');                 

disp('                                                                  ');    
disp('                            THE POSTERIOR MODE IS:                ');
disp('                                                                  ');
disp('   BETA      GAMMA      LAMBDA        Pi*        ZETA_p          NU    RHO_phi      RHO_lambda    RHO_Z  SIGMA_Phi  SIGMA_lambda  SIGMA_Z     SIGMA_R ');
disp(num2str(theta_mode))
disp('                                                                  ');                 
  
  
%=========================================================================
%                          CANDIDATE DENSITY
%=========================================================================

disp('                                                                  ');
disp('            *******STEP 2: HESSIAN AT Mode....    ************'    );
disp('                                                                  ');

Sigma = nhess(objective,param_mode);
Sigma = inv(Sigma);

%=========================================================================
%                            SAVE RESULTS
%=========================================================================

save Matfiles/MH_candidate Sigma param_mode 

path(l);

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;

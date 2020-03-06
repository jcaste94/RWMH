%==========================================================================
%                       DSGE MODEL ESTIMATION 
%                   Metropolis-Hastings Algorithm
%                  (Small NK model in the textbook)
%
% Author: Minsu Chang        minsuc@sas.upenn.edu
% Last modified: 2/24/2016
%
% Edited: Juan Castellanos Silván 
% Date: 2/27/2020
%==========================================================================


%=========================================================================
%                              HOUSEKEEPING
%=========================================================================
clc
clear 
close all
delete *.asv

rng(123); % for replicability of the results

global sample % sample needs to be changed in ``Candidate.m"

tic

l = path;

path('Mfiles',path);
path('Optimization Routines',path);
path('LRE',path);
path('Matfiles',path);

%=========================================================================
%                            EXPERIMENTS
%=========================================================================

sample = 1;         % 1 => Old sample 1983:I - 2002:IV
                    % 2 => New sample 1999:IV - 2019:III

dist = 2;           % 1 => Normal distiribution
                    % 2 => t-Student

                    
disp('                                                                  ');
disp('    BAYESIAN ESTIMATION OF DSGE MODEL: METROPOLIS-HASTINGS        ');
disp('                                                                  ');

%=========================================================================
%                  METROPOLIS-HASTINGS ALGORITHM 
% (Report the Acceptance Rate and Recursive Averages Every 500 draws) 
%=========================================================================

load MH_candidate

Nsim          = input('How Many Posterior Draws?:  ');
disp('                                                                  ');
disp('                                                                  ');

c             = 0.2;
c0            = 0.2;
nu            = 10; % degrees of freedom t-Student
Nburn         = int32(0.50*Nsim)+2;
Thetasim      = zeros(Nsim,length(mode));

% Initialize by taking a draw from a distribution centered at mode
go_on = 0;
while go_on == 0
    if dist == 1
        Thetac = mvnrnd(mode',c0*Sigma);
    else
        Thetac = mvtrnd(c0*Sigma,nu) + mode;
    end
   go_on = (Thetac(8)<=1)*(Thetac(9)<=1)*(Thetac(10)<=1)*(Thetac(2)<=1); % bounds
end
Thetasim(1,:) = Thetac;

accept        = 0;
obj           = dsgeliki(Thetasim(1,:)) + prior(Thetasim(1,:));
counter       = 0;
logposterior  = obj*ones(Nsim,1);
vAcceptanceRate = zeros(Nsim,1);

for i=1:Nsim
    
    if dist == 1
        Thetac = mvnrnd(Thetasim(i,:),c*Sigma);
    else
        Thetac = mvtrnd(c*Sigma,nu) + Thetasim(i,:);
    end
    
    CheckBounds = (Thetac(8)<=1)*(Thetac(9)<=1)*(Thetac(10)<=1)*(Thetac(2)<=1);  % bounds
    
    if CheckBounds == 1 
    
       prioc = prior(Thetac);    
       likic = dsgeliki(Thetac);
       objc  = prioc+likic;       
       
       if objc == -Inf
       
          Thetasim(i+1,:) = Thetasim(i,:);
          logposterior(i+1) = obj;
          
       else % objc > -Inf

          alpha = min(1,exp(objc-obj));
          u = rand(1);

          if u<=alpha
             Thetasim(i+1,:)   = Thetac;
             accept            = accept+1;
             obj               = objc;
             logposterior(i+1) = objc;
          else
             Thetasim(i+1,:)   = Thetasim(i,:);
             logposterior(i+1) = obj;
          end
          
       end % if objc == -Inf
       
    else % CheckBounds NE 1
  
       Thetasim(i+1,:) = Thetasim(i,:);
       logposterior(i+1) = obj;
       
    end  % if CheckBounds == 1

    acceptancerate     = accept/i;
    vAcceptanceRate(i) = acceptancerate;
    counter            = counter + 1;

    if counter==500
       disp('                                                                  ');
       disp(['                               DRAW NUMBER:', num2str(i)]         );
       disp('                                                                  ');
       disp('                                                                  ');    
       disp(['                           ACCEPTANCE RATE:', num2str(acceptancerate)]);
       disp('                                                                  ');
       disp('                                                                  ');    
       disp('                            RECURSIVE AVERAGES                    ');
       disp('                                                                  ');
       disp('   Tau       Kappa      Psi1       Psi2        rA        piA       gammaQ       rho_R       rho_g       rho_z       sigma_R       sigma_g       sigma_z   ');
       disp(num2str(mean(Thetasim(1:i,:))));  
       disp('                                                                  ');
       counter = 0;
    end % if counter==500
    
end %for i=1:Nsim

Thetasim    = Thetasim(Nburn:end,:);

logposterior= logposterior(Nburn:end);

save Matfiles/mhdraws Thetasim logposterior   % Save posterior draws

[Nsim,Npam] = size(Thetasim);

[yy, yy05, yy95]=moment(Thetasim(:,1:13));


%% description
    
%=========================================================================
%                TABLE 1: MEAN AND 5TH-95TH PERCENTILES 
%=========================================================================
if dist == 1
    
    sum_vec = [yy' yy05' yy95'];
    vartype     = {'\tau','\kappa','\psi_1','\psi_2','r^{(A)}',...
                   '\pi^{(A)}','\gamma^{(Q)}',...
                   '\rho_{r}','\rho_{g}', '\rho_{z}', ...
                   '\sigma_{r}','\sigma_{g}', '\sigma_{z}'};

    disp('=========================================================================');
    disp(' Variable Name                       Mean         5%        95%         ');
    disp('=========================================================================');
    for hh=1:length(vartype);
        fprintf('%-30s %10.4f %10.4f %10.4f\n',vartype{hh},sum_vec(hh,1),...
            sum_vec(hh,2),sum_vec(hh,3));    
    end
    disp('========================================================================='); 

    % ----------------------
    % Export table to LaTeX
    % ----------------------
    parameters = {'$\tau$';'$\kappa$'; '$\psi_{1}$';'$\psi_{2}$';'$r^{(A)}$';...
        '$\pi^{(A)}$';'$\gamma^{(Q)}$';'$\rho_{r}$';'$\rho_{g}$'; '$\rho_{z}$';...
        '$\sigma_{r}$'; '$\sigma_{g}$'; '$\sigma_{z}$'};


    T = table(yy', yy05', yy95');
    T.Properties.RowNames = parameters;
    T.Properties.VariableNames{'Var1'} = '\textbf{Mean}';
    T.Properties.VariableNames{'Var2'} = '\textbf{5th perc}';
    T.Properties.VariableNames{'Var3'} = '\textbf{95th perc}';

    path = '/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS4/LaTeX/';
    filename = strcat('tPosteriorEstimates', num2str(sample), num2str(dist),'.tex');
    table2latex(T, strcat(path,filename));

end

%=========================================================================
%                  FIGURE 1: RECURSIVE AVERAGES 
%=========================================================================

pnames = strvcat('\tau','\kappa', '\psi_{1}','\psi_{2}','r^{(A)}',...
    '\pi^{(A)}','\gamma^{(Q)}','\rho_{r}','\rho_{g}', '\rho_{z}', '\sigma_{r}', '\sigma_{g}', '\sigma_{z}');

figure('Position',[20,20,900,600],'Name',...
    'Recursive Averages','Color','w')

rmean = zeros(Nsim,Npam);

for i=1:Nsim
    rmean(i,:) = mean(Thetasim(1:i,:),1);
end

for i=1:(Npam-1)
    
subplot((Npam-1)/3,3,i), plot(rmean(:,i),'LineStyle','-','Color','b',...
        'LineWidth',2.5), grid on, hold on
title(pnames(i,:),'FontSize',12,'FontWeight','bold');    
end


%=========================================================================
%                  FIGURE 2: POSTERIOR MARGINAL DENSITIES 
%=========================================================================

pnames = strvcat('\tau','\kappa', '\psi_{1}','\psi_{2}','r^{(A)}',...
    '\pi^{(A)}','\gamma^{(Q)}','\rho_{r}','\rho_{g}', '\rho_{z}', '\sigma_{r}', '\sigma_{g}', '\sigma_{z}');

figure('Position',[20,20,900,600],'Name',...
    'Posterior Marginal Densities','Color','w')


for i=1:(Npam-1)
    xmin = min(Thetasim(:,i));
    xmax = max(Thetasim(:,i));
    grid = linspace(xmin,xmax,100);
    u    = (1+0.4)*max(ksdensity(Thetasim(:,i)));
subplot((Npam-1)/3,3,i), plot(grid,ksdensity(Thetasim(:,i)),'LineStyle','-','Color','b',...
        'LineWidth',2.5), grid on, hold on
plot([mean(Thetasim(:,i)) mean(Thetasim(:,i))], [0 u],'LineStyle',':',...
    'Color','black','LineWidth',2.5 ), grid on, hold off
axis([xmin xmax 0 u]);
title(pnames(i,:),'FontSize',12,'FontWeight','bold');    
end


x = 29.7;                  % A4 paper size
y = 21.0;                  % A4 paper size
xMargin = 1;               % left/right margins from page borders
yMargin = 1;               % bottom/top margins from page borders
xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[x y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

savepath = '/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS4/LaTeX/';
filename = strcat('pPosteriorMarginalDensities', num2str(sample),num2str(dist),'.pdf');
saveas(gcf, strcat(savepath,filename));


disp('                                                                  ');
disp(['                     ELAPSED TIME:   ', num2str(toc)]             );

elapsedtime=toc;



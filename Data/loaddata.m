%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%                    Importing the data from FRED 
%
%
% Author: Juan Castellanos Silván
% Date: 02/27/2020
%==========================================================================

% Housekeeping
clear all; 
close all; 
clc;

YY = zeros(80,3);

%% Per capita real output growth

% ----------------------
% (1) Real GDP (levels)
% ----------------------
load GDP1.txt
GDP = GDP1(:,2);
clear GDP1


% ----------------
% (2) Population 
% ----------------
load CNP16OV.txt
POP = CNP16OV(:,2);
clear CNP16OV


% ---------------
% GDP per capita
% ---------------
for t = 2:length(GDP)  
    YY(t-1,1) = 100*(log(GDP(t)/POP(t)) - log(GDP(t-1)/POP(t-1)));
end

%% Annualized inflation

load CPIAUCSL.txt
CPI = CPIAUCSL(:,2);
clear CPIAUCSL

for t = 2:length(CPI)  
    YY(t-1,2) = 400*(log(CPI(t)/CPI(t-1)));
end

%% Federal Funds Rate

load FEDFUNDS.txt
FFR = FEDFUNDS(:,2);
clear FEDFUNDS

YY(:,3) = FFR;


%% Write txt file 

% Create a table with the data and variable names
T = table(YY(:,1), YY(:,2), YY(:,3), 'VariableNames', {'GDP', 'CPI', 'FFR'});

% Write data to text file
savepath = '/Users/Castesil/Documents/GitHub/Econ 722 - Schorfheide/PS4/DSGE_Estimation_RWMH/';
name = 'us_update.txt';
writetable(T, strcat(savepath,name),'delimiter','tab')



%==========================================================================
%                       DSGE MODEL ESTIMATION:  
%                Plot the data in the old and new sample
%
%
% Author: Juan Castellanos Silván
% Date: 02/27/2020
%==========================================================================

% Housekeeping
clear all; 
close all; 
clc;

%==========================================================================
%               OLD SAMPLE: 1983:I - 2002:IV
%==========================================================================

% Import 
load us.txt

% Figure
figure(1)
title_names = {'Quarterly Output Growth', 'Quarterly Inflation', 'Federal Fund Rate'};
time = 1983:0.25:2002.75;

for i = 1:size(us,2)
    
  subplot(3,1,i)
  plot(time, us(:,i), 'Color', 'k', 'LineWidth',2)
  title(title_names(i))
  xlim([1983, 2003])
  hold on
  grid on
  
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

path = '/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS4/LaTeX/';
filename ='pDataSample1.pdf';
saveas(gcf, strcat(path,filename));




%==========================================================================
%              NEW SAMPLE: 1999:III-2019:IV
%==========================================================================

% Import 
us_update = load('us_update.txt');

% Figure
figure(2)
title_names = {'Quarterly Output Growth', 'Quarterly Inflation', 'Federal Fund Rate'};
time = 1999.75:0.25:2019.5;

for i = 1:size(us,2)
    
  subplot(3,1,i)
  plot(time, us_update(:,i), 'Color', 'k', 'LineWidth',2)
  title(title_names(i))
  xlim([1999, 2020])
  hold on
  grid on
  
end

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[x y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

path = '/Users/Castesil/Documents/EUI/Year II - PENN/Spring 2020/Econometrics IV/PS/PS4/LaTeX/';
filename ='pDataSample2.pdf';
saveas(gcf, strcat(path,filename));

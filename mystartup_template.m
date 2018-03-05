% Personal start script (stored in ~/.matlab/)
home = getenv('HOME');
repos = getenv('REPOS');   % /home/admin/share/global_bashrc

% this allows you to run matlab from outside the class directory and still see/use the required functions
addpath(strcat(home,'/seismo'));

% ADD Aster library (GEOS 627)
%addpath('/usr/local/matlab_toolboxes/aster/cd_5.2/Lib/');

% ADD GEOTOOLS
%addpath([repos '/GEOTOOLS/matlab_util']);
%startup_geotools

% avoid figures splitting two monitors and going off the top of the screen
% [x0 y0 wid height]
%set(0,'defaultFigurePosition',[50 50 85*8 110*8]);
%set(0,'defaultFigurePosition',[100 300 600 600]);
set(0,'defaultFigurePosition',[100 50 600 600]);

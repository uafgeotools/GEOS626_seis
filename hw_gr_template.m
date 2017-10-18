%
% hw_gr_template.m
% Applied Seismology, GEOS 626, Jan 2014
% 
% Homework template:
% Gutenberg-Richter frequency-magnitude relation
%
% calls seis2GR.m
%

clc
clear
close all

%---------------------------
% Problem 1-1

% load GCMT catalog
% M:     6 x n set of moment tensors
% M0:    1 x n set of scalar seismic moments (derived from moment tensors)
% Mw:    n x 1 set of moment magnitudes (derived from moment tensors)
% dep:   n x 1 set of depths
% eid:   n x 1 set of event IDs
% lat:   n x 1 set of latitudes
% lon:   n x 1 set of longitudes
% otime: n x 1 set of origin times, in seconds (relative to Matlab t=0)
idir = './data/';       % you might need to change the path to data
load([idir 'cmtall_sub']);
whos

% duration of catalog, in years
tran_yr = (max(otime) - min(otime))/365.25;
disp(sprintf('GCMT catalog duration is %.3f years',tran_yr));

% plot the catalog, colored by depth
% note: use wrapTo360 to center the map on longitude=180
[~,isort] = sort(dep,'descend');    % plot deepest events on top
msize = 4^2;                        % marker size
figure;
scatter(wrapTo360(lon(isort)),lat(isort),msize,dep(isort),'filled');
axis([0 360 -90 90]);
xlabel('Longitude'); ylabel('Latitutde');
title(sprintf('GCMT catalog (%i events), colored by depth (km)',length(dep)));
caxis([0 600]); colorbar

% plot deepest events only
%dcut = 600;
%figure;
%plot(wrapTo360(lon(dep > dcut)),lat(dep > dcut),'.');
%axis([0 360 -90 90]);
%xlabel('Longitude'); ylabel('Latitutde');

%---------------------------
% Problem 1-2

dmag = 0.1;     % magnitude bin width
% note: seis2GR.m is NOT a built-in matlab command
[Ncum,N,Medges] = seis2GR(Mw,dmag);

%---------------------------

if 1==1
    % example to show one way to plot multiple items with log-scaled axes
    % (1) generate fake data and a best-fitting line
    % (2) plot lines
    n = 50;
    x = linspace(4,10,n)';
    b = 1;
    a = 9;
    e = 0.1;
    err = e*(-1+2*randn(n,1));
    Nd = 10.^(a-b*x + err);
    N  = 10.^(a-b*x);
    N2  = 10.^(a-1-b*x);

    figure;
    h1 = semilogy(x,Nd,'bV','markersize',8,'markerfacecolor','w');
    hold on;    % after this you do not need semilogy -- it is automatic
    h2 = plot(x,N,'k-');
    h3 = plot(x,N2,'r-');
    legend([h1 h2 h3],'fake data','cumulative fit','incremental fit');
    set(gca,'xtick',4:9);
    set(gca,'ytick',10.^[-1:5],'yticklabel',{'0.1','1','10','100','1000','10000','1000000'});
    xlabel('Moment magnitude, Mw'); ylabel('Number of fake earthquakes');
    title('dummy plot to show one way of plotting multiple things with log-y axes');

    % set iprint=1 if you want to produce figures as postscript files
    % note: type 'ps2pdf plot_example.ps' from the command line to make a pdf
    iprint = 1;
    pdir = './';    % directory where you want the figure to go
    if iprint==1, print(gcf,'-dpsc',[pdir 'plot_example']); end
    
    % an even simpler version (note: no 'hold on' is needed):
    %figure; semilogy(x,Nd,'bV',x,N,'k-',x,N2,'r-');
    %legend('fake data','cumulative fit','incremental fit');
end

%==========================================================================

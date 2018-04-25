%
% tomo_template.m
% GEOS 626, Applied Seismology, Carl Tape
% 
% Template script for homework on seismic tomography.
%

clear
close all
format short
format compact

ax1 = [-121 -114 31 37];        % lon-lat plotting dimensions

%==========================================================================
% LOAD DATA

% load sources
[slon,slat,sind] = textread('events_lonlat.dat','%f%f%f','headerlines',1);
nsrc = length(slat);

% load receivers
[rlon,rlat,rind] = textread('recs_lonlat.dat','%f%f%f','headerlines',1);
nrec = length(rlat);

% load spline centers
q = 8;
[qlon,qlat] = textread('con_lonlat_q08.dat','%f%f','headerlines',0);
nspline = length(qlat);

%==========================================================================
% lon-lat gridpoints for plotting

numx = 100;
[lonplot,latplot,numy,a,b] = gridvec(ax1(1),ax1(2),numx,ax1(3),ax1(4));
nplot = length(lonplot);

% Compute design matrix for expanding a function in terms of splines;
% this is needed to view the tomographic models that we generate at the end.
B = zeros(nplot,nspline);
for ii=1:nspline
    ff = Bkspline(qlon(ii),qlat(ii),q,lonplot,latplot);
    B(:,ii) = ff(:);
end

% EXAMPLE of plotting a function on the 'plotting grid' above.
% the example model vector is is a 1 with all other entries 0
mex = zeros(nspline,1);
mex(96) = 1;
cex = B*mex;       % dimension nplot x 1
% two options
%figure; scatter(lonplot,latplot,4^2,cex,'filled');
figure; pcolor(reshape(lonplot,numy,numx),reshape(latplot,numy,numx),reshape(cex,numy,numx)); shading interp;

%==========================================================================
% INVERSE PROBLEM



%==========================================================================
% PLOTTING THE SOLUTION(S)

% values from GMT 'seis' color palette (type 'colormap(seis)')
seis = [170	0	0;  206	0	0;  243	0	0;  255	24	0;  255	60	0;  255	97	0;
        255	133	0; 255	170	0;  255	206	0;  255	243	0;  255	255	0; 255	255	0;
        231	255	4;  161	255	17;  90	255	30;  51	249	64;  13	242	99;  0	194	152;
        0	125	214;  0	68	248;  0	34	226];
seis = cpt2cmap(seis);



%==========================================================================

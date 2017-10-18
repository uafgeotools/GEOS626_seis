%
% modes_PREMobs.m
%
% Spheroidal and toroidal modes observations listed in PREM (1980).
%
% Preparatory script for sumatra_modes_template.m
%

clear, close all, clc

ddir = './data/';

% spheroidal modes
dfile = [ddir 'prem_Smodes.txt'];
[nobsS,~,lobsS,TS,TstdS] = textread(dfile,'%f%s%f%f%f','headerlines',6);
nS = length(nobsS);
fmhzS = 1./TS*1e3;
isort = [1:nS];
%[~,isort] = sort(fmhzS);   % try this
disp('spheroidal mode observations (measured from seismograms):');
for jj=1:nS
   ii = isort(jj);
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s, f = %6.3f mHz',...
       nobsS(ii),lobsS(ii),TS(ii),TstdS(ii),fmhzS(ii) )); 
end 

% toroidal modes
dfile = [ddir 'prem_Tmodes.txt'];
[nobsT,~,lobsT,TT,TstdT] = textread(dfile,'%f%s%f%f%f','headerlines',6);
nT = length(nobsT);
fmhzT = 1./TT*1e3;
isort = [1:nT];
%[~,isort] = sort(fmhzT);   % try this
disp('toroidal mode observations (measured from seismograms):');
for jj=1:nT
   ii = isort(jj);
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s, f = %6.3f mHz',...
       nobsT(ii),lobsT(ii),TT(ii),TstdT(ii),fmhzT(ii) )); 
end

% plot
% note: 0.2-1.0 mHz is the range of Park et al. (2005), Figure 1
% note: 2S1 had not been observed in 1980
xmin = -0.5; xmax = 10.5; df=0.1; fsize=14;

% toroidal modes (same as in the modes HW)
figure; hold on;
plot([xmin xmax],0.2*[1 1],'k--');
plot([xmin xmax],1.0*[1 1],'k--');
scatter(lobsT,fmhzT,12^2,nobsT,'filled');
colorbar;
xlabel('degree, l'); ylabel('frequency, mHz');
title({'toroidal modes for PREM, colored by n',...
    'open circles = spheroidal modes for PREM'});
axis([xmin xmax 0 4]); caxis([0 3]); grid on;
for ii=1:length(nobsT)
    text(lobsT(ii),fmhzT(ii)+df,sprintf('%iT%i',nobsT(ii),lobsT(ii)),'fontsize',fsize);
end
scatter(lobsS,fmhzS,12^2,nobsS,'ko');
fontsize(14);

% spheroidal modes
figure; hold on;
plot([xmin xmax],0.2*[1 1],'k--');
plot([xmin xmax],1.0*[1 1],'k--');
scatter(lobsS,fmhzS,12^2,nobsS,'filled');
colorbar;
xlabel('degree, l'); ylabel('frequency, mHz');
title({'spheroidal modes for PREM, colored by n',...
    'open circles = toroidal modes for PREM'});
axis([xmin xmax 0 4]); caxis([0 3]); grid on;
for ii=1:length(nobsS)
    text(lobsS(ii),fmhzS(ii)+df,sprintf('%iS%i',nobsS(ii),lobsS(ii)),'fontsize',fsize);
end
scatter(lobsT,fmhzT,12^2,nobsT,'ko');
fontsize(14);

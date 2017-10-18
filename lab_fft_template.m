%
% lab_fft.m
% Carl Tape, Applied Seismology (GEOS 626)
%
% Spectral analysis of the motion of the sun due to planetary orbits.
%
% calls fftvec.m, markT.m
%

clear
close all
format compact

%---------------------------------
% USER PARAMETERS

% switch between using long time series and short time series
longdata = 0;

%---------------------------------

ddir = './data/';

% synthetically derived sample points of the motion of our sun due to planet rotations
load([ddir 'solarshort.dat']);
ts = solarshort(:,1);
xs = solarshort(:,2);
ys = solarshort(:,3);
%cs = [xs + i*ys];
cs = complex(xs,ys);

load([ddir 'solarlong.dat']);
tL = solarlong(:,1);
xL = solarlong(:,2);
yL = solarlong(:,3);
%cL = [xL + i*yL];
cL = complex(xL,yL);
numt = length(tL);

figure; nr=3; nc=2;
slim = 1.5e-3*[-1 1];
ax1 = [slim slim];
ax2 = [0 max(ts) slim];
ax3 = [0 max(tL) slim];

subplot(nr,nc,1); plot(xs,ys,'r'); axis(ax1);
xlabel('x distance, AU'); ylabel('y distance, AU'); axis equal
subplot(nr,nc,2); plot(xL, yL,'b'); axis(ax1);
xlabel('x distance, AU'); ylabel('y distance, AU'); axis equal
subplot(3,1,2); plot(tL, xL,'b',ts,xs,'r--'); axis(ax3);
xlabel('time, years'); ylabel('x distance, AU'); 
subplot(3,1,3); plot(tL, yL,'b',ts,ys,'r--'); axis(ax3);
xlabel('time, years'); ylabel('y distance, AU'); 
%orient tall, wysiwyg

if longdata==0, t = ts; c = cs; else t = tL; c = cL; end

% sidereal orbit periods -- IN YEARS
% http://ssd.jpl.nasa.gov/?planet_phys_par
periods = [0.2408467 0.61519726 1.0000174 1.8808476 11.862615 ...
            29.447498 84.016846 164.79132 247.92065];
np = length(periods);
freqs = 1 ./ periods;
omegs = 2*pi*freqs;
stvec = {'Mer','Ven','Ear','Mar','Jup','Sat','Ura','Nep','Plu'};
disp('ALL TIME UNITS IN YEARS');
disp('          period   frequency  angular frequency');
for ii=1:np
    disp(sprintf('%s %12.6f %11.6f %11.6f',stvec{ii},periods(ii),freqs(ii),omegs(ii)));
end
%          period     frequency  angular frequency
% Mer       0.240847    4.152019   26.087903
% Ven       0.615197    1.625495   10.213286
% Ear       1.000017    0.999983    6.283076
% Mar       1.880848    0.531675    3.340614
% Jup      11.862615    0.084298    0.529663
% Sat      29.447498    0.033959    0.213369
% Ura      84.016846    0.011902    0.074785
% Nep     164.791320    0.006068    0.038128
% Plu     247.920650    0.004034    0.025344

%==========================================================================
% FFT ANALYSIS HERE



%==========================================================================

%
% lab_dispersion.m
% Carl Tape, Applied Seismology (GEOS 626)
% 
% Example of computing group dispersion and phase dispersion using two
% seismograms (real data).
%
% calls bandpass.m, markt.m, fftvec.m
%

clear
close all
clc
format compact

deg = 180/pi;
fsize = 8;

%--------------------------------------------------------------------------
% INPUT PARAMETERS

% target periods for measurements
Ttarvec = [20 30 40 50]';
ftarvec = 1./Ttarvec;
numtar = length(ftarvec);

% distance between PAS and NEE (km)
delx = 331;

% axes limits for dispersion plots of speed (km/s) vs period (s)
ax1 = [18 52 2.8 4.6];
    
% we are told that the phase velocity at each period must fall within these ranges
cran = [3.1 3.9; 3.0 4.3; 3.3 4.5; 3.3 4.5];

%--------------------------------------------------------------------------
% LOAD AND PLOT SEISMOGRAMS

% load data files
file1 = 'pas.dat';
file2 = 'nee.dat';
ddir = './data/';
load([ddir file1]);
load([ddir file2]);

ti = nee(:,1);
dt = ti(2) - ti(1);
ynee = nee(:,2);
ypas = pas(:,2);
% for FFT an even number is extremely helpful, so chop off the first point
ti = ti(2:end);
ynee = ynee(2:end);
ypas = ypas(2:end);
nt = length(ti);

% use hilbert transform to compute envelope
ypasen = abs(hilbert(ypas));
yneeen = abs(hilbert(ynee));

figure; nr=2; nc=1;
xran = [ti(1) ti(end)];

subplot(nr,nc,1); hold on;
plot(ti,ypas,'b');
% plot(ti,ypasen,'k--',ti,-ypasen,'k--');    % envelope
xlabel('Time (s)'); ylabel('Amplitude'); title('Pasadena, LHZ');
xlim(xran);

subplot(nr,nc,2); hold on;
plot(ti,ynee,'r');
% plot(ti,yneeen,'k--',ti,-yneeen,'k--');    % envelope
xlim(xran);
xlabel('Time (s)'); ylabel('Amplitude'); title('Needles, LHZ');
%print(gcf,'-depsc',[pdir 'PAS_NEE_seis']);

%--------------------------------------------------------------------------
% COMPUTE FOURIER TRANSFORM, THEN PLOT AMPLITUDE SPECTRUM

% frequency vector
npt = nt;
f = fftvec(ti);     % note: negative frequencies

Hp = fft(ypas);  %
Ap = abs(Hp);    % =sqrt(H.*conj(H)), where P=H.*conj(H) is the power spectral density
Hn = fft(ynee);
An = abs(Hn);

% explore these to see various details with the matlab FFT
if 0==1
    y = ypas; H = Hp; A = Ap;
    imax = npt/2+1;
    ip = 1:imax;
    
    % note: the first entry contains the INTEGRATED content of the signal
    disp('check the first entry of the FFT:');
    sum(ypas), mean(ypas)*npt, H(1)

    % check the difference between abs and power
    % if z = a + bi, then abs(z) = sqrt(z z*)
    norm( A.*A - H.*conj(H) ) / norm( A )

    % compare IFFT[FFT[y(t)]] with y(t)
    figure; plot(ti,y,'b',ti,ifft(H),'r--');

    % check the ordering of the complex entries of H
    figure; plot(real(H(2:imax-1)) - real(H(npt:-1:imax+1)),'.');
    
    break
end

%-------------
% plot the spectrum for PAS and NEE

% positive and negative frequencies
figure; hold on; 
plot(f,Ap,'b'); plot(f,An,'r');
legend('PAS','NEE');
xlabel('frequency (Hz)'); ylabel('spectral amplitude');
%print(gcf,'-depsc',[pdir 'PAS_NEE_spec_fneg']);

% positive frequencies only
ipos = find(f > 0);
figure; hold on; 
plot(f(ipos),Ap(ipos),'b'); plot(f(ipos),An(ipos),'r');
legend('PAS','NEE');
xlabel('frequency (Hz)'); ylabel('spectral amplitude');
%print(gcf,'-depsc',[pdir 'PAS_NEE_spec']);

%==========================================================================
% CODE HERE FOR GROUP SPEED (use bandpass.m)



%==========================================================================
% CODE HERE FOR HARMONICS (do not use bandpass.m)

% fourier transform of Pasadena and Needles seismograms
Hp = fft(ypas);
Hn = fft(ynee);

% dimensions of variables
whos ti ypas ynee f Hp Hn

tlims = [2580 2780];
for ii=1:numtar
    % target frequency for harmonic
    ftar = ftarvec(ii);
    
    % initialize harmonics
    Hp2 = complex(zeros(npt,1),0);
    Hn2 = complex(zeros(npt,1),0);
    
    % get ftar and -ftar from the frequency vector
    % (this will avoid having matlab tell you that it will ignore the
    % complex conjugate parts when using ifft)
    [~,itemp] = sort(abs(abs(f)-ftar),'ascend');
    itemp = itemp(1:2)
    f(itemp)
    
    % CODE HERE FOR HARMONICS
    %Hp2 
    %Hn2 

end

%==========================================================================
% CODE HERE FOR PHASE SPEED



%==========================================================================

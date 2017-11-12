%
% lab_response.m
% 
% This script shows some conventions for instrument response files
% associated with GISMO/Antelope and rdseed/sac.
% 
% calls numerous GISMO functions and also read_pzfile.m
%
% The example waveform is from CAN (Canberra, Australia) for the
% 2004 Mw 9.X Sumatra-Andaman earthquake.
%
% Carl Tape, 2012-03-31
%

clear, close all, clc

deg = 180/pi;
spdy = 86400;   % seconds per day

% print figures to files
iprint = 0;
pdir = './';

% test file is the 2004 Sumatra earthquake recorded at CAN.G,
% featured in Figure 1 of Park et al. 2005 (Science)
ddir = '/home/admin/databases/SUMATRA/data/';
tdir = [ddir 'CAN_test/'];
dbname = [ddir 'sumatra'];
station = 'CAN';
netwk = 'G';        % note: the network is not added to 'waveform' (get(w,'network'))
channel = 'LHZ';
stlab = [station '_' netwk '_' channel];
pzfile = [tdir 'SAC_PZs_G_CAN_LHZ__1989.153.00.00.00.0000_2006.344.02.60.60.99999'];

% waveform time interval
% note: start time needed to access response file
startTime0 = datenum(2004,12,25,12,58,50);
endTime0 = datenum(2005,1,4,12,58,50);
% pick times to get the entire waveform that is in the database
startTime = startTime0 - 1;
endTime = endTime0 + 1;

dlabs = {'m to counts','m/s to counts','m/s^2 to counts'};

%==========================================================================

% specify frequencies
fmin = 1e-4;
fmax = 1e2;
numf = 100;
f = logspace(log10(fmin),log10(fmax),numf)';
omega = 2*pi*f;
% note: Although f=0 is an allowable frequency for FFT, you will crash the
%       functions below if you have f=0 included. In our homework problems,
%       the frequency vector will be provided when the fft operation is
%       applied.

aran = 10.^[-20 1.5];

% default option: get instrument response from antelope database
% FIGURE 1
res0 = response_get_from_db(station,channel,startTime,f,dbname);
response_plot(res0,[fmin fmax]); ylim(aran);
title('response_get_from_db.m','interpreter','none');
if iprint==1, print(gcf,'-depsc',sprintf('%sCAN_response_fig1',pdir)); end

if 0==1
    % res0 is a matlab structure
    % note: the amplitude of abs(res0.values) is NORMALIZED to 1 (by calib)
    res0
    res0.calib
    res0.units
    Ivcheck = res0.values / (res0.calib*1e-9);  % to match results below
    figure;
    subplot(2,1,1); semilogx(res0.frequencies,angle(Ivcheck));
    subplot(2,1,2); loglog(res0.frequencies,abs(Ivcheck));
    break
end

% Here we create a response object by directly specifying the response file.
% For one file we manually removed the FIR filters in order to show that it
% would then match the response obtained from the converted-to-velocity
% sac pole-zero file.
% FIGURES 2 and 3
for kk=1:2
    if kk==1
        rfile0 = 'STRECKEISEN_STS1.5';
    else
        rfile0 = 'STRECKEISEN_STS1.5_noFIR';
    end
    rfile = [tdir rfile0];
    respObject = dbresponse(rfile);
    response.values = eval_response(respObject,omega);
    response.frequencies = f;
    response_plot(response,[fmin fmax]);
    ylim(aran);
    title(rfile0,'interpreter','none');
    if iprint==1, print(gcf,'-depsc',sprintf('%sCAN_response_fig%i',pdir,kk+1)); end
end

% compare CAN response in antelope database (no FIR) with PZs from sac file
% FIGURE 4 (should match Figure 3)
ideriv = 1;     % response to VELOCITY
listfile = true;
[p,z,c,A0,k] = read_pzfile(pzfile,ideriv,listfile);
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = A0;    % A0 (not c) needed to match normalization
res = response_get_from_polezero(f,polezero);
response_plot(res,[fmin fmax]); ylim(aran);
title(['sac pole-zero file: ' dlabs{ideriv+1}]);
if iprint==1, print(gcf,'-depsc',sprintf('%sCAN_response_fig4',pdir)); end

% relationship between PZ constants (sac) and calib (antelope)
% note: c = A0*k <==> CONSTANT = A0*SENSITIVITY
res0.calib*1e-9, 1/k

%-----------------------
% compare displacement, velocity, and acceleration spectra
% Note that all response plots are normalized such that the
% VELOCITY RESPONSE = 1 at some calibration period.
% FIGURE 5

% first: use read_pzfile.m (add one pole=0 per differentiation)
xf = 5;  % figure index
figure(xf); nr=3; nc=2;
for kk=1:3
    ideriv = kk-1;
    [p,z,c,A0,k] = read_pzfile(pzfile,ideriv);   % c = A0*k
    polezero.poles = p;
    polezero.zeros = z;
    polezero.normalization = c;    % note c, not A0
    res = response_get_from_polezero(f,polezero);
    % complex instrument response to either displacement, velocity, or acceleration
    Ix = res.values;
    % phase response
    subplot(nr,nc,2*kk-1); semilogx(f,angle(Ix)*deg); axis([fmin fmax -180 180]);
    xlabel('frequency, Hz'); ylabel('phase, deg');
    % amplitude response
    subplot(nr,nc,2*kk); loglog(f,abs(Ix)); xlim([fmin fmax]);
    xlabel('frequency, Hz'); ylabel('amplitude');
    [maxI,imax] = max(abs(Ix));
    title({sprintf('sac pole-zero file (%s)',dlabs{kk})...
        sprintf('max = %.3e at %.3e Hz',maxI,f(imax)) });
end

% THIS IS ONLY A CHECK
% second: check this operation by dividing by i*omega for each operation
%         (here we plot in dashed red lines on Figure 10 to show an exact match)
% note: FFT[ f'(t) ] = i*w*FFT[ f(t) ]
[p,z,c,A0,k] = read_pzfile(pzfile,0);   % displacement
polezero.poles = p;
polezero.zeros = z;
polezero.normalization = c;
res = response_get_from_polezero(f,polezero);
Id = res.values;        % displacement
figure(xf); subplot(nr,nc,1); hold on; plot(f,angle(Id)*deg,'r--'); grid on;
figure(xf); subplot(nr,nc,2); hold on; plot(f,abs(Id),'r--'); axis(10.^[-4 2 -5 12]); grid on;
Iv = Id./(1i*omega);    % velocity
figure(xf); subplot(nr,nc,3); hold on; plot(f,angle(Iv)*deg,'r--'); grid on;
figure(xf); subplot(nr,nc,4); hold on; plot(f,abs(Iv),'r--'); axis(10.^[-4 2 -5 12]); grid on;
Ia = Id./(-omega.^2);   % acceleration
figure(xf); subplot(nr,nc,5); hold on; plot(f,angle(Ia)*deg,'r--'); grid on;
figure(xf); subplot(nr,nc,6); hold on; plot(f,abs(Ia),'r--'); axis(10.^[-4 2 -5 12]); grid on;
if iprint==1, orient tall; print(gcf,'-depsc',sprintf('%sCAN_response_fig%i',pdir,xf)); end

%==========================================================================
% EXAMPLES from Alaska

if 0==1
    fmin = 1e-4;
    fmax = 1e3;
    numf = 100;
    f = logspace(log10(fmin),log10(fmax),numf)';

    dbname = '/aerun/sum/params/Stations/master_stations';
    sta = 'F3TN'; chan = 'HHZ';
    res0 = response_get_from_db(sta,chan,datenum(2015,1,1),f,dbname);
    response_plot(res0,[fmin fmax]);
    subplot(2,1,1); title(sprintf('%s : %s %s',dbname,sta,chan),'interpreter','none');
    sta = 'CRP'; chan = 'EHZ';
    res0 = response_get_from_db(sta,chan,datenum(2015,1,1),f,dbname);
    response_plot(res0,[fmin fmax]);
    subplot(2,1,1); title(sprintf('%s : %s %s',dbname,sta,chan),'interpreter','none');
end

%==========================================================================


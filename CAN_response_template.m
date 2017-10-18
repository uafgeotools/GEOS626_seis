%
% CAN_response_template.m
%
% Template script for deconvolving the instrument response in the frequency
% domain.
%
% The example waveform is from CAN (Canberra, Australia) for the
% 2004 Mw 9.X Sumatra-Andaman earthquake.
%
% Be sure to run lab_response.m (and see lab_response.pdf) before proceeding.
%
% Carl Tape, 03/31/2012
%

clear
close all
clc

deg = 180/pi;
spdy = 86400;   % seconds per day

%----------------------------------

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

%==========================================================================
    
% load waveform and remove calibration, so that the raw waveform is IN COUNTS
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,channel,netwk,'');
w = waveform(ds,scnl,startTime,endTime);
wraw = remove_calib(w);
figure; plot(wraw); axis tight
%fontsize(14); print(gcf,'-depsc','./CAN_response_seis');

% get some info about the seismogram
[tstart,tdur,sps] = getm(w,'start','duration','freq');
dt = 1/sps;
disp(sprintf('duration is %.3f days = %.2f hours = %.2f min = %.3e s',...
    tdur,tdur*24,tdur*24*60,tdur*24*60*60));

% example of getting an absolute time from the seismogram
tpick = 3*1e5;      % based on the plot
datestr(tstart + tpick/spdy,31)

% for FFT (specifically for modes): demean and taper
w = demean(w); 
wd = get(w,'data');
nd = length(wd);
taper = tukeywin(nd,1);     % matlab taper function
wd = wd.*taper;
w = set(w,'DATA',wd);
sps = get(w,'freq');        % samples per second
dt = 1/sps;                 % time step

% compute FFT -- this can take 5-10 minutes
% note 1: this will attach the frequency version to the waveform object
% note 2: the file fftmat.m will be saved in your local directory,
%         so be sure to run Matlab from the same directory
fname = 'fftcan';
if ~exist([fname '.mat'],'file')
    tic, w = wf_fft.compute(w); toc
    save(fname,'w');
else
    load(fname);
end
f    = get(w,'fft_freq');
wAmp = get(w,'fft_amp');    % amplitude of C(w)
wPhs = get(w,'fft_phase');  % phase of C(w)
C    = wAmp.*exp(1i*wPhs);  % fourier transform, C(w)

%--------------------------------------------------------------------------
% HOMEWORK/LAB EXERCISE STARTS HERE



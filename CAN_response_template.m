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
% Carl Tape, 2012-03-31
%

clear, close all, clc

deg = 180/pi;
spdy = 86400;   % seconds per day

% load G.CAN.LHZ for Sumatra-Andaman earthquake
station = 'CAN';
channel = 'LHZ';
w = sumatra_read_seis(station,channel);

% example of getting an absolute time from the seismogram
tpick = 3*1e5;      % based on the plot
tstart = get(w,'start');
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

% compute FFT
% NOTES:
% + this can take ~5 minutes, depending on how long the time series is
% + this will attach the frequency version to the waveform object
% + the file fftmat.m will be saved in your local directory,
%         so be sure to run Matlab from the same directory
% + when you copy-and-paste this block, be sure to CHANGE fname (e.g., fftcan_noise)
fname = './fftcan';
if ~exist([fname '.mat'],'file')
    tic, w = wf_fft(w); toc
    save(fname,'w');
else
    load(fname);
end

% compute FFT -- this will attach the frequency version to the waveform object
f    = get(w,'fft_freq');
wAmp = get(w,'fft_amp');    % amplitude of C(w)
wPhs = get(w,'fft_phase');  % phase of C(w)
C    = wAmp.*exp(1i*wPhs);  % fourier transform, C(w)

%--------------------------------------------------------------------------
% HOMEWORK/LAB EXERCISE STARTS HERE



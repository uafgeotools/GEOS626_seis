%
% CAN_response_template.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% Template script for deconvolving the instrument response in the frequency
% domain.
%
% The example waveform is from CAN (Canberra, Australia) for the
% 2004 Mw 9.X Sumatra-Andaman earthquake.
%
% Be sure to run lab_response.m (and see lab_response.pdf) before proceeding.
%

clear, close all, clc

deg = 180/pi;
spdy = 86400;   % seconds per day

% load G.CAN.LHZ for Sumatra-Andaman earthquake
station = 'CAN';
channel = 'LHZ';
w = sumatra_read_seis(station,channel);     % this will also plot the waveform

sps = get(w,'freq');        % samples per second
dt = 1/sps;                 % time step

% example of getting an absolute time from the seismogram
tpick = 3*1e5;      % based on the plot
tstart = get(w,'start');
datestr(tstart + tpick/spdy,31)

% compute FFT
% READE THESE NOTES:
% + this can take ~5 minutes, depending on how long the time series is
% + this will attach the frequency version to the waveform object
% + the file fftmat.m will be saved in your local directory,
%   so be sure to run Matlab from the same directory
% + when you copy-and-paste this block, be sure to CHANGE fname (e.g., fftcan_noise)
fname = './fftcan';
if ~exist([fname '.mat'],'file')
    
    % for FFT (specifically for modes): demean and taper
    w = demean(w); 
    wd = get(w,'data');
    nd = length(wd);
    taper = tukeywin(nd,1);     % matlab taper function
    wd = wd.*taper;
    w = set(w,'DATA',wd);

    tic
    w = wf_fft(w);  % KEY COMMAND
    toc
    
    f    = get(w,'fft_freq');   % frequencies
    wAmp = get(w,'fft_amp');    % amplitude of C(w)
    wPhs = get(w,'fft_phase');  % phase of C(w)
    C    = wAmp.*exp(1i*wPhs);  % fourier transform, C(w)
    
    save(fname,'f','C');
else
    load(fname);
end

whos f C

%--------------------------------------------------------------------------
% HOMEWORK/LAB EXERCISE STARTS HERE



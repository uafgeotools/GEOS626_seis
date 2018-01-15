%
% CAN_bp_template.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
% 
% Template script for bandpass filtering the vertical seismogram.
%

clear, close all, clc

spdy = 86400;       % seconds per day

% load G.CAN.LHZ for Sumatra-Andaman earthquake
station = 'CAN';
channel = 'BHZ';
w = sumatra_read_seis(station,channel);

% specify bandpass filter
% YOU WILL NEED TO CHANGE THE LIMITS HERE
T1 = 100;   % minimum period
T2 = 400;   % maximum period
f1 = 1/T2;
f2 = 1/T1;
npoles = 2;
f = filterobject('B',[f1 f2],npoles);

w = fillgaps(w,'meanAll');
w = demean(w); 
w = detrend(w);
%w = taper(w,0.05);
w0 = w;               % save a copy for later
w = filtfilt(f,w);
figure; plot(w);

% YOU WILL NEED TO APPLY hilbert AND smooth, FOLLOWING NI ET AL. (2005), FIGURE 1
% help waveform/hilbert
% help waveform/smooth
% help smooth
% hint: try something like w = smooth(w,NPT,'moving')
%       where NPT is number of points in the sliding window



%==========================================================================

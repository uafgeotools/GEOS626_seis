%
% CAN_bp_template.m
% 
% Template script for bandpass filtering the vertical seismogram.
%

clear, close all, clc

spdy = 86400;       % seconds per day

% test file is the 2004 Sumatra earthquake recorded at CAN.G, featured in
% Figure 1 of Park et al. 2005 (Science)
ddir = '/home/admin/databases/SUMATRA/data/';
dbname = [ddir 'sumatra'];
station = 'CAN';
channel = 'BHZ';

% waveform time interval
% pick times to get the entire waveform
% NOTE: all times are in days
startTime0 = datenum(2004,12,25,12,58,50);
endTime0 = datenum(2005,1,4,12,58,50);
startTime = startTime0 - 1;
endTime = endTime0 + 1;

% load waveform
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,channel,'','');
w = waveform(ds,scnl,startTime,endTime);
w = remove_calib(w);
datestr(get(w,'start'),31)
datestr(get(w,'end'),31)
get(w,'duration')
figure; plot(w);

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

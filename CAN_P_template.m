%
% CAN_P_template.m
% 
% Template script for analyzing the Sumatra P wave at CAN.
% See also CAN_response_template.m
%

clear
%close all
clc

spdy = 86400;   % seconds per day

% test file is the 2004 Sumatra earthquake recorded at CAN.G, featured in
% Figure 1 of Park et al. 2005 (Science)
ddir = '/home/admin/databases/SUMATRA/data/';
dbname = [ddir 'sumatra'];
station = 'CAN';
channel = {'BHZ','BHE','BHN'};

% waveform time interval
startTime0 = datenum(2004,12,25,12,58,50);
endTime0 = datenum(2005,1,4,12,58,50);
startTime = startTime0 + 12/24;
endTime = startTime + 5000/spdy;

% load waveforms
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,channel,'','');
w = waveform(ds,scnl,startTime,endTime);
w = remove_calib(w);

% plot all three components
figure; nr=3; nc=1;
for kk=1:length(w), subplot(nr,nc,kk); plot(w(kk)); end

% now examine the Z component only
w = w(1);
datestr(get(w,'start'),31)
datestr(get(w,'end'),31)
get(w,'duration')

% example of time picks in seconds
tpicks = [3000 4000];       % CHANGE THESE TO PICK 'THE' P WAVE
tcuts = get(w,'start') + tpicks/spdy;
datestr(tcuts)
figure; plot(w);

% extract a subset time series
% NOTE: this will reset the plotted time axis to start at zero, but the
% absolute time will be correct, as indicated in the title start time;
% you can also plot with the command plot(w,'xunit','date') to show the date
w = extract(w,'TIME',tcuts(1),tcuts(2));
figure; plot(w);

% DO SPECTRAL ANALYSIS HERE


%==========================================================================


%
% CAN_P_template.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
% 
% Template script for analyzing the Sumatra P wave at CAN.
% See also CAN_response_template.m
%

clear, close all, clc

spdy = 86400;   % seconds per day

% load G.CAN for Sumatra-Andaman earthquake
station = 'CAN';
channel = {'BHE','BHN','BHZ'};
w = sumatra_read_seis(station,channel);

% plot all three components
figure; nr=3; nc=1;
for kk=1:length(w), subplot(nr,nc,kk); plot(w(kk)); end

% now examine the Z component only
w = w(3);
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


%
% CAN_noise_template.m
% 
% Template script for analyzing the noise time series at CAN.
% See also CAN_response_template.m
%

clear, close all, clc

spdy = 86400;   % seconds per day

% load the 10-day time series at CAN from BEFORE the Sumatra earthquake
tdir = '/home/admin/databases/SUMATRA/data/sac_LH_noise/G_LH/';
filename = [tdir '2004.350.00.41.42.8698.G.CAN..LHZ.D.SAC'];
w = loadsac(waveform,filename);
figure; plot(w);
tstart = get(w,'start');
tend   = get(w,'end');

% example of getting an absolute time from the seismogram
tpick = 3*1e5;  % CHANGE THIS TO IDENTIFY THE EARTHQUAKE
datestr(tstart + tpick/spdy,31)

% comment this for next step
break

% time picks in seconds
tpicks = [1 8]*1e5;    % CHANGE THESE TO PICK A WINDOW WITH NO VISIBLE EARTHQUAKES
tcuts = get(w,'start') + tpicks/spdy;
datestr(tcuts)

% extract a subset time series
% NOTE: this will reset the plotted time axis to start at zero, but the
% absolute time will be correct, as indicated in the title start time;
% you can also plot with the command plot(w,'xunit','date') to show the date
w = extract(w,'TIME',tcuts(1),tcuts(2));
figure; plot(w);

% SPECTRAL ANALYSIS AND DECONVOLUTION


%==========================================================================


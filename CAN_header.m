%
% CAN_header.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% This script reads in a pre-saved seismic waveform the shows the 2004
% Sumatra-Andaman earthquake recorded at Canberra, Australia (G.CAN).
% The seismogram was featured in Figure 1 of Park et al. 2005 (Science).
%

%==========================================================================
% origin time for Sumatra-Andaman earthquake

% www.globalcmt.org produces a few different output formats for earthquakes
% in the GCMT catalog
%
% A: Standard
% This format lists the centroid time.
%
% 122604A OFF W COAST OF NORTHERN
% 
%   Date: 2004/12/26   Centroid Time:  1: 1: 9.0 GMT
%   Lat=   3.09  Lon=  94.26
%   Depth= 28.6   Half duration=95.0
%   Centroid time minus hypocenter time: 139.0
%   Moment Tensor: Expo=29  1.040 -0.427 -0.610 2.980 -2.400 0.426 
%   Mw = 9.0    mb = 8.9    Ms = 8.9   Scalar Moment = 3.95e+29
%   Fault plane:  strike=329    dip=8   slip=110
%   Fault plane:  strike=129    dip=83   slip=87
%
% B: CMTSOLUTION
% This format lists the PDE origin time at the top and the time shift to
%    centroid time = PDE origin time + tshift
%
%  PDE 2004 12 26  0 58 50.00   3.3000   95.7800  10.0 8.9 8.9 OFF W COAST OF NORTHERN                 
% event name:     122604A        
% time shift:    138.9600
% half duration:  95.0000
% latitude:        3.0900
% longitude:      94.2600
% depth:          28.6100
% Mrr:       1.040000e+29
% Mtt:      -4.270000e+28
% Mpp:      -6.100000e+28
% Mrt:       2.980000e+29
% Mrp:      -2.400000e+29
% Mtp:       4.260000e+28
%

otimePDE    = datenum(2004,12,26,0,58,50);
tshiftCMT_s = 139;
otimeCMT    = otimePDE + tshiftCMT_s/86400;

disp(sprintf('Sumatra-Andaman PDE origin time    : %s',datestr(otimePDE,31)));
disp(sprintf('Sumatra-Andaman GCMT centroid time : %s',datestr(otimeCMT,31)));

%==========================================================================

% if multiple channels are listed, just used the first one here
if length(channel) > 1, chanx = channel(1); else chanx = char(channel); end

% waveform time interval
% pick times to get the entire waveform that is saved in the database
% note: this will add NaN to before and after, if no waveform is available
if strcmp(chanx,'BHZ')
    startTime = otimePDE - 0.14/24;
    endTime = startTime + 2.05/24;
else
    startTime = otimePDE - 0.5;
    endTime = startTime + 10;
end

% test file is the 2004 Sumatra earthquake recorded at CAN.G,
% featured in Figure 1 of Park et al. 2005 (Science)
ddir = '/home/admin/databases/SUMATRA/data/';
dbname = [ddir 'sumatra'];
station = 'CAN';
netwk = 'G';        % note: the network is not added to 'waveform' (get(w,'network'))
 
% load waveform and remove calibration, so that the raw waveform is IN COUNTS
ds = datasource('antelope',dbname); 
scnl = scnlobject(station,chanx,netwk,'');
w = waveform(ds,scnl,startTime-0.001,endTime);   % -0.001 is a kludge
d = get(w,'data');
disp(sprintf('number of NaN points is %i',sum(isnan(d))));

% check for NaN points
% (for some reason, there is always a aingle NaN point at the start)
w = extract(w,'time',startTime,endTime);
d = get(w,'data');
disp(sprintf('number of NaN points is %i',sum(isnan(d))));

w = remove_calib(w);
figure; plot(w); axis tight
%fontsize(14); print(gcf,'-depsc','./CAN_response_seis');

% get some info about the seismogram
[tstart,tdur,sps] = getm(w,'start','duration','freq');
dt = 1/sps;
disp(sprintf('duration is %.3f days = %.2f hours = %.2f min = %.3e s',...
    tdur,tdur*24,tdur*24*60,tdur*24*60*60));

%==========================================================================

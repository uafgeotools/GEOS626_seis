function [w,otimePDE,otimeCMT] = sumatra_read_seis(station,channel)
%SUMATRA_READ_SEIS read seismogram from database of Sumatra-Andaman waveforms
%
% This script reads in a pre-saved seismic waveform of the 2004
% Sumatra-Andaman earthquake.
%
% EXAMPLES:
%    w = sumatra_read_seis('CAN','BHZ');
%    w = sumatra_read_seis('CAN','LHE');
%    w = sumatra_read_seis('QSPA','BHZ_00');
%    w = sumatra_read_seis('CAN',{'BHE','BHN','BHZ'});
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

startTime = otimePDE - 1;
endTime   = otimePDE + 10;

% test file is the 2004 Sumatra earthquake recorded at CAN.G,
% featured in Figure 1 of Park et al. 2005 (Science)
ddir = '/home/admin/databases/SUMATRA/data/';
dbname = [ddir 'sumatra'];
 
% load waveform and remove calibration, so that the raw waveform is IN COUNTS
netwk = '';     % network is not a critical input for the Sumatra database
ds = datasource('antelope',dbname); 
% WARNING: If you ask for channel = {'BHZ','BHE','BHN'},
%          the waveforms will be returned in the order of E-N-Z, not Z-E-N.
%          It is safest to provide the input as channel = {'BHE','BHN','BHZ'},
scnl = scnlobject(station,channel,netwk,'');
w = waveform(ds,scnl,startTime,endTime);        % KEY COMMAND

nw = length(w);

for ii=1:nw 
    check_nan(w(ii));    

    % remove NaN points
    % (for some reason, there is always a single NaN point at the start)
    d = get(w(ii),'data');
    igood = find(~isnan(d)==1);
    w(ii) = extract(w(ii),'index',igood(1),igood(end));
    check_nan(w(ii));

    % remove calibration, so that the raw waveform is in COUNTS
    w(ii) = remove_calib(w(ii));
    figure; plot(w(ii));
    %axis tight    % will help reveal NaN in time series
    %fontsize(14); print(gcf,'-depsc','./CAN_response_seis');

    % get some info about the seismogram
    tdur = getm(w(ii),'duration');
    disp(sprintf('duration is %.3f days = %.2f hours = %.2f min = %.3e s',...
        tdur,tdur*24,tdur*24*60,tdur*24*60*60));
end

% if multiple channels are selected, this will ensure that all have the
% same time window
if nw > 1
    [tstart,tend] = getm(w,'start','end');
    tstart0 = max(tstart);
    tend0 = min(tend);
    w = extract(w,'time',tstart0,tend0);
    tdur = getm(w(1),'duration');
    disp(sprintf('duration is %.3f days = %.2f hours = %.2f min = %.3e s',...
        tdur,tdur*24,tdur*24*60,tdur*24*60*60));
end

%==========================================================================

function check_nan(w)

d = get(w,'data');
disp(sprintf('number of NaN points is %i',sum(isnan(d))));

%==========================================================================

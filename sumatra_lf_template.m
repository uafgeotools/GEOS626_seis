%
% sumatra_lf_template.m
% 
% Template script for analyzing direct arrival waveforms from Sumatra:
%    channel LHZ, duration up to 10 days
%
% This assumes that you have added the path to the GEOTOOLS directories.
%

clear, close all, clc

spdy = 86400;   % seconds per day

% extract the full database of BHZ waveforms
otimePDE = datenum(2004,12,26,0,58,50);
originTime = otimePDE;
tshift = 0.5*3600;      % in seconds
duration_s = 4.0*3600;  % in seconds
startTime = originTime - tshift/spdy;
endTime   = originTime + duration_s/spdy;
dur_dy = endTime-startTime;
fprintf('startTime is %s\n',datestr(startTime,31));
fprintf('total length of time requested: %.2f s (= %.2f min = %.2f hours)\n',...
    dur_dy*spdy,dur_dy*3600,dur_dy*24);

% the advantage of using getwaveform.m is that it will add all kinds of
% headers to the waveform objects, such as station azimuth, distance, etc
idatabase = 7;
channel = {'LHZ'};
iint = 0;
iprocess = 1;           % calibration applied (nm/s)
cutoff = []; samplerate = []; stasub = []; sacdir = [];
% source parameters from GCMT catalog
elat = 3.09;
elon = 94.26;
edep_km = 28.6;
mag = 8.9974;
eid = 'M122604A';
% KEY COMMAND to getwaveform.m
w = getwaveform(idatabase,startTime,endTime,channel,iint,...
    iprocess,cutoff,samplerate,stasub,sacdir,originTime,elat,elon,edep_km,mag,eid);
nw = length(w);

disp('here is a list of the waveforms you have:');
for ii=1:nw
   disp(sprintf('%3i %7s %3s %6s %10s',ii,get(w(ii),'channel'),get(w(ii),'KNETWK'),...
       get(w(ii),'station'),get(w(ii),'units')));
end

% save a copy to avoid rerunning
w0 = w;

%------------------------------------------
% RERUN THE BLOCK BELOW WHEN SELECTING YOUR STATIONS

% pick a subset of waveforms
% note: alternatively you can use wkeep.m
%ipick = [1:nw];                % default
ipick = [22 30 31 120:124];     % CHANGE THIS     
w = w0(ipick);

% PLOTTING PARAMETERS FOR plotw_rs.m (CHANGE THESE AS NEEDED)
rssort = 2;      % =1 by azimuth, =2 by distance
iabs = 0;
tmark = [];
pmax = 40;
iintp = 0;
inorm = 1;
tlims = [];
nfac = 1;
azstart = [];
iunit = 2;
imap = 0;

% plot record section
T1 = 100;
T2 = 1000;
plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

% plot map
[sta,rlat,rlon,elat,elon] = getm(w,'station','STLA','STLO','EVLA','EVLO');
plot_event_station(elat(1),elon(1),rlat,rlon,sta);
%------------------------------------------

% START YOUR ANALYSIS HERE


%==========================================================================
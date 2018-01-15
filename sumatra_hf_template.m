%
% sumatra_hf_template.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
% 
% Template script for analyzing the direct arrival waveforms from Sumatra:
%     channel BHZ, duration of up to 2 hours
%
% This assumes that you have added the path to the GEOTOOLS directories.
%

clear, close all, clc

spdy = 86400;   % seconds per day

if 0==1
    % quick plot example
    startTime = 7.323055408564815e+05;
    endTime = 7.323175408564815e+05;
    ds = datasource('antelope','/home/admin/databases/SUMATRA/data/sumatra');
    scnl = scnlobject('CAN','LHZ','G','');
    w = waveform(ds,scnl,startTime,endTime);
    figure; plot(w,'xunit','h');
end

% extract the full database of BHZ waveforms
otimePDE   = datenum(2004,12,26,0,58,50);
originTime = otimePDE;
startTime  = originTime - 1/24;
endTime    = originTime + 4/24;

% the advantage of using getwaveform.m is that it will add all kinds of
% headers to the waveform objects, such as station azimuth, distance, etc
idatabase = 7;
channel = {'BHZ'};
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

% pick a subset of waveforms
% note: alternatively you can use wkeep.m
%ipick = [1:nw];                    % default
ipick = [21 22 107 67 68 165];     % USER: CHANGE THIS     
w = w0(ipick);

% PLOTTING PARAMETERS FOR plotw_rs.m (USER: CHANGE THESE AS NEEDED)
rssort = 2;      % =1 by azimuth, =2 by distance
iabs = 0;
tshift = [];
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
T1 = [];
T2 = [];
plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

% plot global map
[sta,rlat,rlon,elat,elon] = getm(w,'station','STLA','STLO','EVLA','EVLO');
plot_event_station(elat(1),elon(1),rlat,rlon,sta);

% SOME EXAMPLES OF USING THE PLOTTING COMMANDS

% example of cutting a record
% note: you can alternatively use wcut.m
w([4 6]) = [];  % cut a bad record (WAKE) and a repeated record (KDAK)
plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

% OPTION A: example of applying a relative time shift to each seismogram
% Note: This is in the order of listed stations (NOT as ordered in the record section).
% Note: The DT are w.r.t. the origin time and are listed on the labels.
get(w,'station')
tshift = [1186 1250 845 1440];
plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

% an easier way (though less accurate) is to estimate the P travel time
% based on the source-station distance
dist_deg = get(w,'GCARC');
if 0==1
    % OPTION B: assume a constant velocity (as suggested in the homework)
    Vest = 11;                      % km/s
    Ptt = deg2km(dist_deg) / Vest;  % very crude estimation for 30 < Delta < 85
    
else
    % OPTION C: use the Jeffreys-Bullen traveltime tables
    % WARNING: our simplified version only works for epicentral distances <100 deg
    Ptt = get_JB_Ptime(edep_km,dist_deg);
    if any(isnan(Ptt))
        disp('WARNING: JB times are NaN, since no direct P for Delta > 100 deg');
        Ptt(isnan(Ptt)) = 1000;  % dummy arrival time (use tauP in the future!)
    end
    for ii=1:length(w), w(ii) = addfield(w(ii),'JBP',Ptt(ii)); end
    % check that you added a new field to each w(ii)
    get(w,'JBP')
end

% now replot
% Here the time shift is relative to the marker time (originTime),
% so the DT in the record section is the predicted P travel time (Ptt).
tshiftmark = Ptt;
tmark = originTime;
plotw_rs(w,rssort,iabs,tshiftmark,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

% example of resetting plotting range
% Note that the amplitude scaling is based on the full-length seismogram,
% not the (subset) time interval that is plotted.
tlims = [-50 700];
plotw_rs(w,rssort,iabs,tshiftmark,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);

% START YOUR ANALYSIS HERE


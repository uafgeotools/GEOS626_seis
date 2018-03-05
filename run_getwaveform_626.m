%
% run_getwaveform_626.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% Load waveforms from the UAF waveform database, then plot record sections.
%
% copied from run_getwaveform_short.m in GEOTOOLS
% Abbreviated version of run_getwaveform.m and getwaveform_input.m
%

%clear              % do not clear variables when replotting record sections 
clc, close all

%==========================================================================
% USER PARAMETERS

iex = 1;                % CHANGE THIS OR ADD YOU OWN EXAMPLE
bgetwaveform = false;    % 
bplotrs = true;         % 

%==========================================================================

% simple example of getting waveforms from a database
if 0==1
    startTime = datenum(2012,4,11,8,41,57);
    endTime   = datenum(2012,4,11,10,21,57);
    ds = datasource('uaf_continuous');
    sta = {'MDM'};
    cha = {'BHZ'};
    scnl = scnlobject(sta,cha,'','');
    w = waveform(ds,scnl,startTime,endTime);
    figure; plot(w);
end

spdy = 86400;

% default parameters for plotting record sections
% note: you can over-ride these within each example below
rssort = 2;      % =1 by azimuth, =2 by distance
iabs = 0;
T1 = [];
T2 = [];
trshift = 0;
tmark = [];
pmax = 50;
iintp = 0;
inorm = 1;
nfac = 1;
azstart = [];
iunit = 1;
tlims = [];     % time limits for plotting
imap = 1;

%==========================================================================

% default parameters for waveform extraction
samplerate = [];
cutoff = [];

switch iex
    case 1
        % Yahtse event recorded within two databases (AEC, Yahtse)
        %idatabase = 5; stasub = [-141.489 -141.377 60.205 60.245];
        idatabase = [1 5]; stasub = [0 200];

        % source parameters (some can be empty)
        originTime = datenum('2010/09/18 14:15:02');
        elat = 60.155496;
        elon = -141.378343;
        edep_km = 0;
        eid = [];
        mag = [];

        chan = {'HHZ','BHZ'};

        duration_s = 70;
        oshift = 20;
        T1 = 0.1;
        T2 = 2;

    case 2
        % Mw 7.5 SE Alaska
        idatabase = 1; stasub = [0 500];
        % source parameters (some can be empty)
        originTime = datenum('2013/01/05 08:58:32.4');
        elat = 55.62;
        elon = -135.13;
        edep_km = 0;
        eid = [];
        mag = [];
        % broadband channels
        chan1 = {'BHZ','BHE','BHN','BH1','BH2'};
        % strong motion channels
        chan2 = {'BNZ','BNE','BNN','BN1','BN2','BLZ','BLE','BLN','BL1','BL2',...
            'HNZ','HNE','HNN','HN1','HN2','HLZ','HLE','HLN','HL1','HL2'};
        % warning: waveforms will have different units (nm/s, nm/s^2)
        chan = [chan1 chan2];
        duration_s = 300;
        oshift = 50;
        
    case 3
        % explosion in Fairbanks
        idatabase = 1;
        % source parameters (some can be empty)
        originTime = datenum('2013/02/03 01:10:31.495');
        %elat = 64.8156; elon = -147.9419;  % original AEC
        %elat = 64.8045; elon = -147.9653;  % reviewed AEC
        elat = 64.80175; elon = -147.98236; % infrasound
        edep_km = 0;
        eid = [];
        mag = [];
        % broadband channels
        chan = {'SHZ','HHZ','BHZ'};
        stasub = [0 200];
        duration_s = stasub(2)/0.30;  % air wave
        oshift = 50;
        % record section
        T1 = 0.2;
        T2 = 1;
        
    case 4
        % very low frequency (VLF) event near Kantishna
        idatabase = 1;
        % source parameters (some can be empty)
        originTime = datenum('2014/01/22 12:14:34');
        elat = 63.463; elon = -150.109;
        edep_km = 38;
        eid = [];
        mag = 1.7;
        % broadband channels
        chan = {'BHZ'};
        stasub = [0 200];
        duration_s = 100; 
        oshift = 0;
        T1 = []; T2 = 2;
        
    case 5
        % landslide near Lituya
        idatabase = 1; stasub = [0 1000];
        % source parameters (some can be empty)
        %originTime = datenum('2010/09/18 14:15:02');
        originTime = datenum('2014/02/16 14:24:00');
        elat = 58.68;
        elon = -137.37;
        edep_km = 0;
        eid = [];
        mag = [];
        chan = {'BHZ'};  % consider HHZ as well
        duration_s = 600;
        oshift = 0;
        T1 = 10; T2 = 40;
        pmax = 40;
        
    case 6
        % Sumatra Mw 8.6 in Alaska and triggered earthquakes
        idatabase = 1;
        % four different events
        originTimeS = datenum('2012/04/11 08:38:37.30');    % Sumatra (CMT-PDE)
        originTimeA = datenum('2012/04/11 09:00:09.71');    % Andreanof (NEIC)
        originTimeN = datenum('2012/04/11 09:21:57.44');    % Nenana (NEIC)
        originTimeI = datenum('2012/04/11 09:40:58.02');    % Iliamna slab (NEIC)
        elatS = 2.24; elonS = 92.78; edepS = 40.03; emagS = 8.6;
        elatA = 51.36; elonA = -176.10; edepA = 20; emagA = 5.5;
        elatN = 64.9222; elonN = -148.9461; edepN = 19.4; emagN = 3.88;
        elatI = 60.10; elonI = -152.83; edepI = 101; emagI = 2.9;
        
        % source
        elat = elatS; elon = elonS; edep_km = edepS; mag = emagS;
        eid = '201204110838A';  % GCMT
        originTime = originTimeS;
        chan = {'BHZ'};
        dmax_deg = 25;
        stasub = [elonN elatN 0 dmax_deg 0 360];

        % parameters for Sumatra waveforms across Alaska
        duration_s = 3600*2.0;
        oshift = 3600*0.25;
        T1 = 1/4; T2 = 1/2;     % P wave + triggered events
        %T1 = 2; T2 = 1000;     % full wavetrain (no triggered events visible)
        iunit = 2;
        % three different triggered events are visible
        tmark = [originTimeA originTimeN originTimeI];
        
    case 7
        % one of the triggered earthquakes
        idatabase = 1;
        % four different events
        originTimeS = datenum('2012/04/11 08:38:37.30');    % Sumatra (CMT-PDE)
        originTimeA = datenum('2012/04/11 09:00:09.71');    % Andreanof (NEIC)
        originTimeN = datenum('2012/04/11 09:21:57.44');    % Nenana (NEIC)
        originTimeI = datenum('2012/04/11 09:40:58.02');    % Iliamna slab (NEIC)
        elatS = 2.24; elonS = 92.78; edepS = 40.03; emagS = 8.6;
        elatA = 51.36; elonA = -176.10; edepA = 20; emagA = 5.5;
        elatN = 64.9222; elonN = -148.9461; edepN = 19.4; emagN = 3.88;
        elatI = 60.10; elonI = -152.83; edepI = 101; emagI = 2.9;
        
        % source: USER CHOOSE ONE
        originTime = originTimeN; elat = elatN; elon = elonN; edep_km = edepN; mag = emagN; stasub = [0  200]; duration_s = 200;  % Nenana
        %originTime = originTimeA; elat = elatA; elon = elonA; edep_km = edepA; mag = emagA; stasub = [0 2000]; duration_s = 600;  % Andreanof
        %originTime = originTimeI; elat = elatI; elon = elonI; edep_km = edepI; mag = emagI; stasub = [0  400]; duration_s = 200;  % Iliamna
        eid = [];
        chan = {'BHZ'};

        % bandpass
        oshift = 10; T1 = 1/4; T2 = 1/2;
        
    case 8
        % TRY YOUR OWN EXAMPLE HERE
        
        
end

% tshift (and trshift) are for plotting record sections only
tshift = oshift + trshift;

%startTime = originTime - max(oshift)/spdy;
startTime = originTime - oshift/spdy;
endTime   = originTime + duration_s/spdy;
dur_dy = endTime-startTime;
fprintf('origin time is %s\n',datestr(originTime,'yyyy-mm-dd HH:MM:SS.FFF'))
fprintf('startTime is %s\n',datestr(startTime,31));
fprintf('total length of time requested: %.2f s (= %.2f min = %.2f hours)\n',...
    dur_dy*spdy,dur_dy*3600,dur_dy*24);

% additional user parameters
%sacdir = './';      % =[] to return waveform object only
sacdir = [];
iint = 0;            % integrate waveforms: =1 for displacement, =0 for velocity
iprocess = 1;        % iprocess = 2 to deconvolve
irs = 1;

if bgetwaveform
    % get waveform object (optional: write sac files to a directory)
    tic
    [w,s,site,sitechan] = getwaveform(idatabase,startTime,endTime,chan,iint,iprocess,cutoff,samplerate,stasub,sacdir,originTime,elat,elon,edep_km,mag,eid);
    toc
    disp(sprintf('%.1f s to execute getwaveform.m from run_getwaveform_short.m',toc));
end

%==========================================================================
% FROM HERE ON OUT, PLOTTING AND CHECKING ONLY

if bplotrs

whos w s site sitechan

% check indexing and units
for ii=1:length(w)
   disp(sprintf('%3i %7s %3s %6s %6s %6s %10s %6i sps',ii,get(w(ii),'channel'),get(w(ii),'KNETWK'),...
       get(w(ii),'station'),char(site(ii).sta),char(sitechan(ii).sta),get(w(ii),'units'),get(w(ii),'freq')));
end

% plot record section
if and(irs==1,~isempty(w))
    % assume all waveforms have the same event ID
    keid = get(w(1),'KEVNM');

    % plot record section
    % note: these can be printed to file by setting bprint_record_section=true in plotw_rs.m
    plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);
end

end  % bplotrs

%--------------------------------------------------------------------------
% ANALYZE INDIVIDUAL RECORDS

break

% subset of stations
stasubset = {'BOOM'};       % use BOOM with example 1
wsubset = wkeep(w,stasubset);

% enter the following settings
% 1. click AutoUpdate
% 2. slide Freq Max to the maximum
% 3. set NFFT to 256
% 4. set OVERLAP to 80%
% 5. set MIN 42
% 6. set MAX 112
% 7. from the tab at top, Plot --> Specgram2
% WHEN YOU ARE DONE OR IF YOU WANT TO MAKE A NEW PLOT, CLOSE THE FIGURE
figure; uispecgram(wsubset(1));     % will only work for a single 

% in case you want to analyze the filtered version that you see in the
% record section, you can return wfilt
[isort,wfilt] = plotw_rs(w,rssort,iabs,tshift,tmark,T1,T2,pmax,iintp,inorm,tlims,nfac,azstart,iunit,imap);
get(w,'station')
get(wfilt,'station')
get(wfilt(isort),'station')
get(wfilt,'station')
% plot kth station that appears in the ORDERED record section
kk = 1;     % order 
% plot that single record
figure; plot(wfilt(isort(kk)));

%==========================================================================

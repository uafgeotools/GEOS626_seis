%
% sumatra_modes_template.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
% 
% Template script for analyzing modes spectra for sumatra
%
% For a view of the 139/169 spectra from 'kept' seismograms, open this file:
%    /home/admin/databases/SUMATRA/data/wfobject/all_sumatra_modes.pdf
%
% This assumes that you have added the path to the GEOTOOLS directories.
%

clear
close all
clc

bload = true;      % USER: CHANGE THIS

% add path in order to access additional files and scripts
datadir = '/home/admin/databases/SUMATRA/data/wfobject/';
addpath(datadir);

% first load all the data
if bload
    % get waveforms from the database (calls getwaveform.m)
    sumatra_modes_fft;
    disp(sprintf('%i/%i bad records that will not be used:',length(scut),nw));
    stas(scut)
    
    % plot the 30 time series that were manually cut from the analysis
    for ii = 1:nw   % loop over all 169 time series
        if any(ii==scut)  
            stag = [stas{ii} '_' chans{ii} '_' nets{ii}];
            stdur = sprintf('duration = %.2f days',get(w(ii),'duration'));
            stit = [stag ', ' stdur];
            disp(sprintf('%i/%i %s',ii,nw,stag));
            figure; hold on; plot(w(ii)); title(stit);
        end
    end
    
    break
end

% load list of files
% note: column 1 is the index into the full set of FFT waveforms [1:169]
%       column 2 is the index into the reduced set of FFT waveforms;
%                this also represents the page number of all_sumatra_modes.pdf
[ind,ind_pdf,sta,chan,net,tag,ikeep] = textread([datadir 'sumatra_modes.txt'],'%f%f%s%s%s%s%f');
disp('ALL AVAILABLE WAVEFORMS');
for ii=1:length(ind)
   disp(sprintf('%3i %3i %7s %7s %4s %4i',ind(ii),ind_pdf(ii),sta{ii},chan{ii},net{ii},ikeep(ii)));
end
% subset of files with pre-computed FFTs
sta  = sta(logical(ikeep));
chan = chan(logical(ikeep));
net  = net(logical(ikeep));
tag  = tag(logical(ikeep));
ind_pdf = ind_pdf(logical(ikeep));
disp('SUBSET OF WAVEFORMS FOR SCIENTIFIC ANALYSIS (all_sumatra_modes.pdf)');
for ii=1:length(ind_pdf)
   disp(sprintf('%3i %7s %7s %4s',ind_pdf(ii),sta{ii},chan{ii},net{ii}));
end

% USER: PICK A SET TO PLOT AND SAVE FOR ANALYSIS
% note: these are the same indices as the page numbers of the composite PDF
%ipick = [1:length(sta)];    % view all 139 not-cut stations
ipick = [1:3 21];           % (random default)

npick = length(ipick);
w(1,npick) = waveform;  % initialize array of waveforms
                        % WARNING: w might contain an empty waveform object
                        %          if you ask for one
pamp = zeros(npick,1);
for ii=1:length(ipick)
    jj = ipick(ii);
    % load PRE-STORED seismogram and frequency spectrum
    % note: these were saved in sumatra_modes_fft.m / sumatra_loop.m
    % note: time series had calibration applied,
    %       but the full instrument response was NOT removed
    fname = strcat('w',tag{jj},'.mat');
    ifile = [datadir 'full_length/' fname];
    load(ifile);
    w(ii) = v;

    f = get(v,'fft_freq');
    A = get(v,'fft_amp');
    figure; plot(f*1e3,A); xlim([0.2 1.0]);
    title(tag{jj},'interpreter','none');
    xlabel('frequency, mHz'); ylabel('amplitude');

    % save amplitude of a particular peak
    f1 = 0.94; f2 = 0.95;   % USER: CHANGE THESE
    pamp(ii) = max(A(and(f*1e3 > f1, f*1e3 < f2)));
end

% diplay the properties of the object
w(4)
% get properties of the object
% here are some examples
% STLO  station longitude
% STLA  station latitude
% GCARC source-station epicentral distance, in degrees
% AZ    source-station azimuthal angle, in degrees measured from north
[slon,slat,dist_deg,az] = getm(w(4),'STLO','STLA','GCARC','AZ');

break

% example of stacking some spectra
f1 = 0.2*1e-3;
f2 = 10*1e-3;
n = 1e5;
[Astack,fplot,A] = w2fstack(w,f1,f2,n);
figure; plot(fplot*1e3,Astack);
xlabel('frequency, mHz');
ylabel('stacked amplitude spectrum');
% note: matrix A contains all spectra
figure; plot(fplot,A);

% START MODES PROBLEM HERE


%==========================================================================

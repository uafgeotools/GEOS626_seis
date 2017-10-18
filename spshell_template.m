%
% spshell_template.m
% Applied Seismology, Geos 626, UAF
%
% This program computes the toroidal modes for a spherical shell with
% uniform rigidity and density.
%
% calls these functions:
%    stress_disp_tor.m
%    surf_stress.m
%
% GLOBAL VARIABLES: l rvec WT rspan mu rho
%
% by Charles Ammon, Pen State, 2000
% modifications by Carl Tape, UAF, 01/2012
%

close all; clear all; clc

% global variables
% WARNING: DO NOT CHANGE THE DIMENSION OF ANY OF THESE VARIABLES
global l rvec WT rspan imod rho mu % omega

%------------------------------------------------
% USER INPUT

% shell spans from core-mantle boundary (b) to surface (a)
earthr = 6371000;       % radius of earth, in meters (a)
cmbr = 3480000;         % radius of core-mantle boundary, in meters (b)
rspan = [cmbr earthr];  % [b a]

% Earth model
imod = 0;               % index for Earth model (see earthfun.m)
                        %   =0 for homogeneous
                        %   =1 for linear rho(r) and mu(r)
                        %   =2 for cubic rho(r) and mu(r)
% WARNING: BECAUSE THESE ARE GLOBAL VARIABLES, THESE VALUES MAY BE
% OVER-RIDDEN BY VALUES DEFINED IN OTHER FUNCTIONS.                        
rho = 4380;             % density
mu  = 5930*5930*rho;    % rigidity (mu = 1.54e11 Pa)

% options for searching solution space
l = 2;                  % degree (l >= 1)
nmax = 8;               % maximum n (default = 8)
                        % nmax+1 is the max number of roots/eigenfunctions/subplots
                        % nmax=0 will return the first root (n=0)
% plotting options                        
iplot_eig_freqs = 1;    % plot eigenfunctions (=1) or not (=0)
                        % =0 will speed up the calculations
iplot_all_freqs = 1;    % plot eigenfunctions for all input frequencies,
                        % including those that do not satisfy the boundary conditions
                        % =1 for initial problem only

% path to the directory containing the data file prem_Tmodes.txt
ddir = './data/';

iprint = 0; % print figures to file (=1) or not (=0)
pdir = './';

%------------------------------------------------

% range of frequencies (note: omega = 2*pi*f), in Hz
fmin = 1/3600;      % initial frequency to start (T = one hour)
df = 0.0002;        % frequency step size (chosen by trial and error)
%fmax = 0.08;        % stopping frequency (somewhat arbitrary)
fmax = 0.003;
fvec = [fmin:df:fmax];
numf = length(fvec);
disp(sprintf('frequency vector ranges from %.3f mHz to %.3f mHz',fmin*1e3,fmax*1e3));
disp(sprintf('num frequency points is %i, df = %.3f mHz',numf,df*1e3));
disp(sprintf('--> period ranges from %.2f min to %.2f min',1/fmin/60,1/fmax/60));

% THIS BLOCK IS FOR INITIAL PLOTTING ONLY
if iplot_all_freqs==1
    for ii=1:numf
        disp(sprintf('%2i/%2i: f = %.3f mHz',ii,numf,fvec(ii)*1e3));
        % update W(r) and T(r), stored within WT
        surf_stress(fvec(ii));
        
        % plotting parameters
        rplot = rvec/1000;
        xmx = 1.1; ymn = rspan(1)/1000; ymx = rspan(2)/1000; dy = 100;

        % displacement for each frequency
        Wplot = WT(:,1)/max(abs(WT(:,1)));
        figure(2); hold on; plot(Wplot,rplot,'b');
        text(Wplot(end),rplot(end)+dy/2,num2str(ii));
        plot([0 0],rspan/1000,'k','linewidth',2);
        xlabel('normalized displacement, W(r)'); ylabel('radius, km');
        axis([-xmx xmx ymn-dy ymx+dy]);

        % stress for each frequency
        Tplot = WT(:,2)/max(abs(WT(:,2)));
        figure(3); hold on; plot(Tplot,rplot,'r');
        text(Tplot(end),rplot(end)+dy/2,num2str(ii));
        plot([0 0],rspan/1000,'k','linewidth',2);
        xlabel('normalized stress, T(r)'); ylabel('radius, km');
        axis([-xmx xmx ymn-dy ymx+dy]);
    end
    
    % print figures for HW
    if iprint==1
        figure(2); print(gcf,'-depsc',[pdir 'modes_Wr']);
        figure(3); print(gcf,'-depsc',[pdir 'modes_Tr']);
    end

    % exit
    break
end

% initial freqeuncy and corresponding surface stress
% NOTE: surf_stress.m calls stress_disp_tor.m, which depends on degree l
ii = 1;
f = fvec(ii);
Tsurf = surf_stress(f);
n = 0;             % counter for n (n=0 is the first root)

disp(sprintf('freq %3i/%3i %10.3e %10.3e %.2f mHz %.1f s %.2f min', ...
        ii, numf, NaN, Tsurf, f*1e3, 1/f, 1/f/60));

% leave gap for T(n=0,l=1), which do not exist
% note: these are useful when looping over degree l
if and(l==1,nmax>0), n = 1; end        % fill the n >= 1 entries
if and(l==1,nmax==0), continue; end     % exit loop early (mode 0T1 does not exist)

% THIS IS THE KEY LOOP OVER FREQUENCIES
froots = NaN*ones(nmax+1,1);
for ii = 2:numf
    % frequency interval over which we check for a root
    oldf = f;
    f = fvec(ii);
    
    % The function surf_stress.m will updated the key variable WT,
    % which contains the radial displacement W(r) in the first column
    % and stress T(r) in the second column.
    Tsurfold = Tsurf;          % surface stress for previous f
    Tsurf = surf_stress(f);    % surface stress for new f
    
    disp(sprintf('freq %3i/%3i %10.3e %10.3e %.2f mHz %.1f s %.2f min', ...
        ii, numf, Tsurfold, Tsurf, f*1e3, 1/f, 1/f/60));

    % Check if the value of the surface-stress has changed sign,
    % which would indicate that we passed at least one root.
    % If we did cross a root, call fzero to refine the root.
    % Then store the root in the vector froots and plot the results.
    if (Tsurfold * Tsurf < 0)
        % note: fzero is a built-in matlab function
        f0 = fzero('surf_stress',[oldf f]);
        froots(n+1) = f0;

        % update eigenfunctions (WT, rvec) for the exact frequency
        surf_stress(f0);
        disp(sprintf('T(a)=0 --> n=%i %.3f mHz l=%i', n,f0*1e3,l));

        % plotting eigenfunctions (displacement and stress as a function of radius)
        if iplot_eig_freqs==1
            xmx = 1.2; ymn = rspan(1)/1000; ymx = rspan(2)/1000;
            rplot = rvec/1000;
            Wplot = WT(:,1)/max(abs(WT(:,1)));
            Tplot = WT(:,2)/max(abs(WT(:,2)));
            figure(1); if nmax==0, subplot(1,1,n+1); else subplot(3,3,n+1); end
            hold on;
            plot(Wplot,rplot,'b');    % W(r), displacement (blue)
            plot(Tplot,rplot,'r');    % T(r), stress (red)
            plot([-xmx xmx],ymn*[1 1],'k',[-xmx xmx],ymx*[1 1],'k',[0 0],[ymn ymx],'k');
            axis([-xmx xmx ymn-300 ymx+300]); %grid on;
            title(sprintf('f = %.2f mHz, T = %.2f min (l = %i)',froots(n+1)*1000,1/f0/60,l));
            text(-1,ymx+100,sprintf('n = %i',n));
            if mod(n-1,3)==0, ylabel('radius (km)'); end
        end
        
        % exit the loop when you reach n=nmax
        if n==nmax, break; end
        % count for the next root
        n = n + 1;
    end
end
fprintf('/// l = %i, nroots = %i (nmax = %i, fmax = %.3f mHz)\n',...
    l,sum(~isnan(froots)),nmax,fmax*1e3);

break

% toroidal mode observations used in PREM
dfile = [ddir 'prem_Tmodes.txt'];
[nobs,~,lobs,Tobs,Tobs_std] = textread(dfile,'%f%s%f%f%f','headerlines',6);
disp('normal mode observations (measured from seismograms):');
for ii=1:length(nobs)
   disp(sprintf('n = %i, l = %2i, T = %8.2f +/- %.2f s',nobs(ii),lobs(ii),Tobs(ii),Tobs_std(ii))); 
end
% it may be useful to store these in a matrix
maxl = max(lobs);
maxn = max(nobs);
Tmat = NaN(maxn+1,maxl);
for ii=1:length(nobs)
   Tmat(nobs(ii)+1,lobs(ii)) = Tobs(ii);
end

if 0==1
    % example code for calculating misfit in arrays with NaN
    % say that A is observations, B is predictions
    A = rand(4,10); A([3 5 21 3]) = NaN;
    B = rand(4,10); B([3 12 14]) = NaN;
    a = A(:);
    b = B(:);
    r = a - b;
    % you cannot sum a vector that has NaN
    sum(r.^2)
    % so just sum the values that are NOT NaN
    dmis = sum(d(~isnan(d)).^2)
end

%==========================================================================

%
% genlsq.m
%
% genlsq = Generalized Least Squares
% Template script for the iterative quasi-Newton method for a 4-parameter
% inversion for epicenter (xs, ys), origin time (ts), and velocity (V).
% The algorithm employs generalized least squares, where by both data
% covariances and model covariances are used.
%
% Background reading: Tarantola book (2005), Ch. 3 and Appendix 6.22
%
% calls forward_epicenter.m, plot_epicenters.m, plot_covsamples.m
%
% Carl Tape, 01-Jan-2012
%

clc
clear
close all
format compact
format short

% add path to Aster library (chi2inv)
%(if you need this and if it is not already in your path)
%addpath('/usr/local/matlab_toolboxes/aster/cd_5.2/Lib/');

%=========================================
% USER INPUT

nsamples = 1000;
irandom_initial_model = 0;      % 0(fixed), 1(random)
irandom_target_model = 0;       % 0(fixed), 1(random)
idata_errors = 2;               % 0(none),  1(random), 2(fixed)
ifig = 1;                       % 0,1

% to print figures to PS files
iprint = 0;
pdir = pwd;

%=========================================

inormalization = 1;
stnsamples = [num2str(nsamples) ' samples'];
stlabS = {'Sd(m^k)','Sm(m^k)','S(m^k) = Sd + Sm'};
   
%---------------------------------------------
% FORWARD PROBLEM    

forward_epicenter;      % KEY COMMAND

% predictions for prior and initial models (not necessary)
dprior   = d(mprior);
dinitial = d(minitial);

if ifig==1
    % plot different histograms of properties of the prior model covariance samples
    figure; nr=2; nc=2;
    for kk=1:nparm
        sigma = sigma_prior(kk);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = cov_samples_m(kk,:);
        subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.4]); grid on;
        title({['mprior samples: Model parameter ' num2str(kk) ' (' mlabs{kk} ')'],
            ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ]});
    end
    if iprint==1, print(gcf,'-depsc',[pdir 'mprior1']); end

    % plot different histograms of properties of the data covariance samples
    figure; nr=4; nc=3;
    for ii=1:ndata
        sigma = sigma_obs(ii);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = cov_samples_d(ii,:);
        subplot(nr,nc,ii); plot_histo(etemp,edges); ylim([0 0.4]);
        title({['Data index ' num2str(ii)],
            ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ]});
    end
    orient tall, wysiwyg
    if iprint==1, print(gcf,'-depsc',[pdir 'CD1']); end
    
    figure; hold on;
    plot(dobs_samples,'.-');
    p1 = plot(dprior,'bo-','linewidth',2,'markersize',10,'markerfacecolor','b','markeredgecolor','w');
    p2 = plot(dinitial,'ko-','linewidth',2,'markersize',10,'markerfacecolor','k','markeredgecolor','w');
    p3 = plot(dtarget,'ro--','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
    p4 = plot(dobs,'ro-','linewidth',2,'markersize',10,'markerfacecolor','r','markeredgecolor','w');
    legend([p1 p2 p3 p4],'g(mprior)','g(minitial)','g(mtarget)','g(mtarget) + errors',...
        'location','northwest');
    %title(' BLACK = d(mprior);  RED DASHED = d(mtarget);  RED = d(mtarget) + errors');
    xlim([0.5 ndata+0.5]); set(gca,'xtick',[1:ndata]);
    xlabel('Data index'); ylabel('Prediction value, g(m)'); grid on;
    if iprint==1, print(gcf,'-depsc',[pdir 'CD2b']); end
end

%---------------------------------------------
% MISFIT FUNCTION : least squares, Tarantola (2005), Eq. 6.251
% (This calls the function d to compute the predictions.)

% data misfit
Sd = @(m,dobs,icobs) ( 0.5* ...
      transpose(d(m)-dobs) * icobs * (d(m)-dobs) );
% model misfit (related to regularization)
Sm = @(m,mprior,icprior) ( 0.5* ...
      transpose(m - mprior) * icprior * (m - mprior) );
% total misfit
S = @(m,dobs,mprior,icobs,icprior) ( Sd(m,dobs,icobs) + Sm(m,mprior,icprior) );

% initial model
%mnew = mprior;     % prior model
mnew = minitial;
Sd_0 = Sd(mnew,dobs,icobs);
Sm_0 = Sm(mnew,mprior,icprior);
S_0  = S(mnew,dobs,mprior,icobs,icprior);
stS0 = sprintf(' S(m0) = %.3f = %.3f(D) + %.3f(M)',S_0,Sd_0,Sm_0);
disp(stS0);

niter = input(' Select the number of iterations (< 10) or 0 to exit: ');
if niter==0, break; end

% initialize arrays
iter_vec = [0:niter]';
Sd_vec = zeros(niter+1,1);
Sm_vec = zeros(niter+1,1);
S_vec = zeros(niter+1,1);

% misfit for initial model
Sd_vec(1) = Sd_0;
Sm_vec(1) = Sm_0;
S_vec(1) = S_0;

for nn = 1:niter
    %///////////////////////////////
    % CODE HERE for quasi-Newton algorithm

    

    % fill misfit function S_vec, Sd_vec, Sm_vec for plotting later
    %Sd_vec(nn+1) = 
    %Sm_vec(nn+1) = 
    %S_vec(nn+1) = 
    
    %///////////////////////////////
end

% misfit function values
disp('summary of misfit function:');
disp(sprintf('%8s%16s%16s%16s','iter','Sd','Sm','S = Sm + Sd'));
for nn = 1:niter+1
    disp(sprintf('%8i%16.10f%16.10f%16.10f',iter_vec(nn),Sd_vec(nn),Sm_vec(nn),S_vec(nn)));
end

if ifig==1
    % plot convergence curve
    ylims = 10.^[-2 2];
    stit = [num2str(niter) ' iterations'];
    figure; hold on;
    plot(iter_vec, log10(Sd_vec),'r.-',iter_vec, log10(Sm_vec),'b.-',iter_vec, log10(S_vec),'k.-',...
        'linewidth',2,'markersize',20);
    legend(stlabS); xlim([-0.5 niter+0.5]); ylim(log10(ylims));
    set(gca,'xtick',[-1:niter+1]); grid on;
    xlabel('k, iteration'); ylabel(' log10[ S(m^k) ], misfit function'); title(stit);
    if iprint==1, print(gcf,'-depsc',[pdir 'converge_' ftag]); end
end

%///////////////////////////////
% COMPUTE THE FOLLOWING
% mpost       posterior model ("final" model)
% dpost       predictions for mpost
% Gpost       partial derivatives matrix at mpost
% cpost0      posterior covariance matrix (use icobs0 and icprior0)
% sigma_post  variances of the posterior covariance matrix
% rho_post    posterior correlation matrix (hint: see Tarantola, Section 3.3)
% CODE HERE




%///////////////////////////////

% a priori model correlations (for comparison)
rho_prior = corrcov(cprior0);

cpost0, rho_post

% posterior data covariance matrix (e.g., Tarantola Eq. 3.44)
cpost0_d = Gpost*cpost0*Gpost';
cpost0_d = (cpost0_d + cpost0_d')/2;    % force to be symmetric
sigma_post_d = sqrt( diag(cpost0_d) );  % IGNORING OFF-DIAGONAL ELEMENTS
rho_post_d = corrcov(cpost0_d);         % posterior correlation matrix
rho_prior_d = corrcov(cobs0);           % prior, for comparison

%format long
disp(sprintf('model summary (%i iterations):',niter));
disp('    prior    initial  posterior   target');
disp([mprior minitial mpost mtarget]);
disp(sprintf('data summary (%i observations):',ndata));
disp('    prior    initial   posterior   target   actual');
disp([dprior dinitial dpost dtarget dobs]);

% Cholesky decomposition to obtain the square-root of cpost0
% NOTE: for large problems, this is not possible due to poor
%       conditioning of cpost0 or the inability to compute cpost0
Lpost = chol(cpost0,'lower');

% samples of the posterior distribution
mpost_samples = zeros(nparm,nsamples);
mcov_samples = zeros(nparm,nsamples);
randn_vecs_m = randn(nparm,nsamples);
mcov_samples  = Lpost * randn_vecs_m;
mpost_samples = repmat(mpost,1,nsamples) + mcov_samples;

% compare the standard deviation with sigma_post
std_samples = std(mpost_samples');

% compare posterior model distribution with prior
% note: format statement allows for vectors (like sigma_prior)
disp('  ');
disp(' Compare model uncertainties : ');
disp(['             model parameter : ' sprintf('%13s', mlabs{:}) ]);
disp(['                       units : ' sprintf('%13s', ulabs{:}) ]);
disp(['                 sigma_prior = ' sprintf('%13.5f', sigma_prior) ]);
disp(['                  sigma_post = ' sprintf('%13.5f', sigma_post) ]);
disp(['   std(' sprintf('%6.0f', nsamples) ' mpost_samples) = ' sprintf('%13.5f', std_samples) ]);
disp(['    sigma_post / sigma_prior = ' sprintf('%13.5f', sigma_post./sigma_prior) ]);
disp('  ');

% compute the predictions associated with the posterior samples,
% then compare std_d_samples with sigma_post_d
d_samples = zeros(ndata,nsamples);
for xx=1:nsamples
    ms = mpost_samples(:,xx);
    d_samples(:,xx) = d(ms);
end
covd_samples = cov(d_samples');
rhod_samples = corrcov(covd_samples);
std_d_samples = sqrt(diag(covd_samples));
%std_d_samples = std(d_samples'); 

disp('  ');
disp(' Compare data uncertainties : ');
disp(sprintf('%16s %10s %10s %10s','prior','post','samples','post/prior'));
call = [sigma_obs(:) sigma_post_d(:) std_d_samples(:) sigma_post_d(:)./sigma_obs(:)];
for ii=1:ndata
   disp(sprintf('%6i%10.4f %10.4f %10.4f %10.4f',ii,call(ii,:)));
end
disp('  ');

if ifig==1
    % display distributions for each model parameter (nparm ROWS of cov_samples_m)
    figure; nr=2; nc=2;
    for kk=1:nparm
        sigma = sigma_post(kk);
        edges = [-4*sigma: sigma/2 : 4*sigma];
        etemp = mcov_samples(kk,:);
        subplot(nr,nc,kk); plot_histo(etemp,edges); ylim([0 0.4]); grid on;
        stl1 = 'mpost samples';
        stl2 = ['Model parameter ' num2str(kk) ' (' mlabs{kk} ')'];
        stl3 = ['mean = ' sprintf('%.5f',mean(etemp)) ...
            '; std = ' sprintf('%.5f',std(etemp)) ];
        if kk==1, title({stl1,stl2,stl3});
        else title({stl2,stl3}); end
    end
    if iprint==1, print(gcf,'-depsc',[pdir 'mpost2_' ftag]); end
    
    % correlation matrices and scatterplots
    %plot_covsamples(mprior_samples,rho_prior,'mprior',[],[],[],mlabs);
    plot_covsamples(mpost_samples,rho_post,'mpost',[],[],[],mlabs);
    plot_covsamples(mprior_samples,rho_prior,'mprior',mpost_samples,rho_post,'mpost',mlabs);
    
    % 'physical view' of estimated posterior data uncertainties
    % note: plot either sigma_post_d (from Cpost_d) or std_d_samples (from d(Cpost_samples))
    figure; hold on;
    plot(mpost_samples(1,:),mpost_samples(2,:),'c.');
    plot(mpost(1),mpost(2),'o','markersize',10,'markerfacecolor','c','markeredgecolor','w');
    %plot(mprior(1),mprior(2),'o','markersize',10,'markerfacecolor','b','markeredgecolor','w');
    %scatter(xrec,yrec,16^2,sigma_post_d,'filled','V'); title('estimated uncertainties for posterior predictions'); 
    scatter(xrec,yrec,16^2,std_d_samples,'filled','V'); title('uncertainties for posterior predictions, computed from samples'); 
    scatter(xrec,yrec,16^2,'kV');
    colormap('hot'); colorbar;
    axis equal; axis(axepi);
    %set(gca,'xtick',[0:20:100],'ytick',[0:20:100]);
    xlabel('X distance (km)'); ylabel('Y distance (km)');
    
    % plot predictions for samples of the posterior
    plot_covsamples(d_samples,rho_post_d,'dpost',[],[],[],[]);
    
    % opts is set in forward_epicenter
    plot_epicenters(mprior_samples,mprior,minitial,mtarget,opts,mpost);
    % plot the cpost0 samples and re-plot the two markers
    plot(mpost_samples(1,:),mpost_samples(2,:),'c.');
    plot(mpost(1),mpost(2),'o','markersize',10,'markerfacecolor','c','markeredgecolor','w');
    plot(mtarget(1),mtarget(2),'o','markersize',10,'markerfacecolor','r','markeredgecolor','w');
    title('samples of prior (blue) and posterior (cyan)');
    if iprint==1, print(gcf,'-depsc',[pdir 'mpost_' ftag '_epi']); end
end

%==========================================================================

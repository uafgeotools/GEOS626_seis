%
% lab_linefit.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% This program introduces the least squares method for the example of
% fitting a line (i.e., a model with two parameters) to a set of scattered
% data.  We show three options for solving the problem, each which gives
% the same result. No regularization is needed (or used) in this example.
%
% See lab_linefit.pdf for instructions to this exercise.
%

clear
close all
%format long
%format compact

%---------------------------

% USER PARAMETERS (CHANGE THESE)
n = 50;                     % number of observations
sigma = 0.3;                % standard deviations of added errors

% TARGET model vector (y-intercept, slope)
mtar = [2.1 -0.5]';
m = length(mtar);           % number of model parameters

%---------------------------

% compute design matrix
% x is a vector of x_i locations where your measurements y_i will be made
xmin = -2;
xmax = 2;
x = linspace(xmin,xmax,n)'; % input x-values
G = [ones(n,1) x];          % n by m design matrix
                         
% display dimensions of these variables
whos

% generate errors
e = sigma * randn(n,1); % normally distributed random numbers

% generate target 'data' with errors added
dtar = G*mtar;
d = dtar + e;

% optional: add one big anomaly
%d(1) = d(1) + 1000*sigma;

% SOLVE: compute least squares solution, estimates, and estimated variance.
% (We show several options for mest, each with the same result.)
mest = G\d
%mest = inv(G'*G)*G'*d;             % note matlab's warning about using inv
%mest = pinv(G)*d;
%mest = flipud(polyfit(x,d,1)')     % specific to this problem only

dest = G*mest;          % estimated predictions
res = d - dest;         % residuals

figure; msize = 10;
stres = [' std(res) = ' sprintf('%.3f', std(res) )];

subplot(2,1,1); hold on;
plot(x,d,'.','markersize',msize);
plot(x,dtar,'c--','linewidth',3);
plot(x,dest,'r--','linewidth',3);
legend('data','mtar','mest','location','southwest');
xlabel(' x'); ylabel(' d');
title({sprintf('Estimated model : m = (%.2f, %.2f)',mest(1),mest(2)), stres})
grid on; axis equal;
axis([min(x) max(x) min(G*mtar)-2*sigma max(G*mtar)+2*sigma]);

subplot(2,2,3);
plot(res,'.','markersize',msize); grid on; ylim([-1 1]);
xlabel(' Observation index'); ylabel(' Residual, d - dest'); title(stres);

subplot(2,2,4);
edges = [-1.05:0.1:1.05]; [Nh,bin] = histc(res,edges);
bar(edges,Nh,'histc'); xlim([min(edges) max(edges)]);
xlabel(' Residual'); ylabel(' Number'); title([' Ntotal = ' num2str(n)]);

%fontsize(11); orient tall, wysiwyg

%error

%---------------------------
% generate a plot showing the residual sum of squares (RSS) as a function of model space

% search range, measured by the distance from the target model
m1_ran = 1;
m2_ran = 1;

% number controlling the number of gridpoints in model space
nx = 100;   % for the misfit function plot
%nx = 10;   % for the gradient plot (fewer arrows)

% generate grid for the model space
m1_vec = linspace(mtar(1)-m1_ran, mtar(1)+m1_ran, nx);
m2_vec = linspace(mtar(2)-m2_ran, mtar(2)+m2_ran, nx);
[M1,M2] = meshgrid(m1_vec,m2_vec);
[a,b] = size(M1);
ng = a*b;                       % number of gridpoints in model space
m1 = reshape(M1,1,ng);
m2 = reshape(M2,1,ng);

% compute misfit function (and gradient)
RSSm = zeros(ng,1);             % initialize misfit function
% INITIALIZE GRADIENT HERE


for kk=1:ng
    mtry = [m1(kk) m2(kk)]';    % a sample from model space
    dtry = G*mtry;              % predictions from the model
    res = d - dtry;             % residuals between data and predictions
    RSS = sum(res.*res);        % residual sum of squares
    RSSm(kk) = RSS;
    
    % COMPUTE GRADIENT HERE [DO NOT USE MATLAB'S gradient FUNCTION]

end

% plot the misfit function
nc = 30;    % number of contours to plot
figure; hold on;
contourf(M1,M2,reshape(RSSm,a,b),nc); shading flat
%pcolor(M1,M2,reshape(RSSm,a,b)); shading flat
%scatter(m1,m2,6^2,RSSm,'filled'); shading flat;
l1 = plot(mtar(1),mtar(2),'ws','markersize',10,'markerfacecolor','k');
l2 = plot(mest(1),mest(2),'wo','markersize',10,'markerfacecolor','r');
legend([l1,l2],'target model','estimated model');
axis equal, axis tight;
caxis([-1e-6 0.5*max(RSSm)]); colorbar
xlabel('m1, y-intercept');
ylabel('m2, slope');
title('Residual sum of squares');

% PLOT GRADIENT HERE WITH quiver COMMAND ('help quiver')


%==========================================================================

%
% lab_newton.m
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% This is a template for implementing an iterative Newton algorithm.
% See Tarantola (2005), Eq. 6.291.
% See lab_newton.pdf for instructions to this lab exercise.
%

close all, clear, clc
format long

% define a misfit function and its derivatives
% note: our 'model vector' is one-dimensional
P = [2 -1 2 -1 2];          % non-quadratic
%P = [0  0 2 -1 2];         % quadratic
%F = @(m) ( polyval(P,m) );
F = @(m) ( P(1)*m.^4 + P(2)*m.^3 + P(3)*m.^2 + P(4)*m + P(5) );
g = @(m) ( 4*P(1)*m.^3 + 3*P(2)*m.^2 + 2*P(3)*m + P(4) );
H = @(m) ( 12*P(1)*m.^2 + 6*P(2)*m + 2*P(3) );
whos

% specify bounds for choosing initial model (and for plotting)
mA = -2;
mB = 2;
mplot = linspace(mA,mB,100);

% COMPUTE MINIMUM OF MISFIT FUNCTION HERE USING BUILT-IN MATLAB FUNCTION


figure; hold on;
plot(mplot,F(mplot));
% PLOT MATLAB MINIMUM HERE
xlabel('model m'); ylabel('misfit function, F(m)');

%---------------------------------

% initial model
m0 = -1.5;
%m0 = mA + (mB-mA)*rand; % try a random starting value

% IMPLEMENT NEWTON ALGORITHM HERE (see class notes "Taylor series and least squares")



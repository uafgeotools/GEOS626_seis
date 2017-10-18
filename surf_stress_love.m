function value = surf_stress_love(k0)
% modified from surf_stress.m in the modes homework

global k rspan rvec WT mbeta mmu omega

% this updates k (global) for stress_disp_love.m
% note: we need surf_stress_love.m to have an input variable in order to
%       apply fzero to find its roots; here we use k0 since we do not want
%       to confuse our global k with the local k0
k = k0; 

% calculate initial conditions at r=0 within the mantle halfspace
mk = omega/mbeta;
nub = sqrt(k^2 - mk^2);
if ~isreal(nub), warning('setting nub=0 (k=%.3e mk=%.3e)',k,mk); nub=0; end
Tbot = mmu*nub;
WT0 = [1.0 Tbot];   % the initial values of [displacement stress]

% integrate the system of differential equations y' = f(t,y)
% this calls our function stress_disp_love.m to
% calculate derivatives (d/dr) of displacement and stress;
% on return the vectors rvec and WT contain the values of radius
% and the displacement and stress eigenfunctions
% note: the dimension of rvec and WT is the number of points needed for
%       the numerical integration -- this will vary

% try this option to use default tolerance on numerical solution
[rvec,WT] = ode45('stress_disp_love',rspan,WT0);

% try this option if you need really accurate eigenfunctions
%opts = odeset('reltol',1e-12,'abstol',1e-12);
%[rvec,WT] = ode45('stress_disp_love',rspan,WT0,opts);

value = WT(end,2);      % stress value at earth's surface (r = rspan(2))

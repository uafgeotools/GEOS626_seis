function dWT = stress_disp_tor(r,WT)
%STRESS_DISP_TOR displacement and stress for a toroidal earth
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% The input values of WT(1) and WT(2) are W(r) and T(r) respectively.
% The returned deriatives are stored in dWT

global omega l imod rho mu % rspan rvec WT

% structural values at radius r: density and rigidity
% note: if imod=0, then the program will use the rho and mu from spshell.m
if imod~=0, [rho,mu] = earthfun(r); end

% displacement (first row of equation 1)
dWT(1,1) = WT(1) / r + WT(2) / mu;

% stress (second row of equation 2)
dWT(2,1) = ((l-1)*(l+2)*mu/(r*r) - rho*omega*omega)*WT(1) - 3*WT(2)/r;

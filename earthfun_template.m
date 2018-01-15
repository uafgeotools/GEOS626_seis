function [rho,mu] = earthfun(r)
%EARTHFUN return a rho and mu value for a specified radius r
% Applied Seismology, GEOS 626, University of Alaska Fairbanks
%
% called by stress_disp_tor.m

global rspan imod

switch imod
    case 1
        % linear model
        % ENTER YOUR CODE HERE
        cmbr = rspan(1);    % b
        earthr = rspan(2);  % a
        
        error('earthfun.m imod=1 not yet implemented (comment this out when you have implemented it)');
        
    case 2
        % cubic model
        % ENTER YOUR CODE HERE
        
        error('earthfun.m imod=2 not yet implemented (comment this out when you have implemented it)');
    
    otherwise
        error('invalid imod (=1,2)');
end

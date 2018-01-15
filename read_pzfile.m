function [p,z,c,A0,k] = read_pzfile(pzfile,ideriv,listfile)
%READ_PZFILE reads sac-convention poles and zeros from a sac pole-zero file
%
% INPUT
%   pzfile  pole-zero file produced by rdseed (v 5.2)
%   ideriv  =0 for disp (default), =1 for vel, =2 for accel
%   iwrite  OPTIONAL: if present, write the PZ file to the command window
%
% There are three essential points about the sac PZ files that are produced
% by the program rdseed.
%    1. The response is always for DISPLACEMENT in units of meters.
%    2. Filters associated with the digitizer or decimation (FIR) are
%          NOT included (though they are in the RESP file).
%    3. The INPUT is meters, the OUTPUT is counts (not the other way).
%
% note: CONSTANT = A0 * SENSITIVITY (c = A0*k)
% note: the antelope 'calib' field is 1/CONSTANT (adjusting for units)
%
% Carl Tape, 2012-03-29
%

if ~exist(pzfile,'file'), error([pzfile ' does not exist']); end
lines = textread(pzfile,'%s','delimiter','\n','whitespace','');

% write file to command window
if nargin==3
    if listfile
        disp('======== read_pzfile.m ========');
        for ii=1:length(lines), disp(lines{ii}); end
        disp('======== read_pzfile.m ========');
    end
end

for ii=1:length(lines)
    ltemp = lines{ii};
    ch1 = ltemp(1);
    % selected header information
    if ch1=='*'
       if strcmp(ltemp(3:4),'A0')
          A0 = str2double(strtrim(ltemp(22:end)));
       end
       if strcmp(ltemp(3:13),'SENSITIVITY')
          k = str2double(strtrim(ltemp(22:34)));
       end
    end
    % zeros
    if ch1=='Z'
        nz = str2double(strtrim(ltemp(6:end)));
        z = zeros(nz,2);
        for kk=1:nz
           %ltemp0 = lines{ii+kk};
           %z(kk,1) = str2double(strtrim(ltemp0(1:14)));
           %z(kk,2) = str2double(strtrim(ltemp0(15:end)));
           z(kk,:) = str2num(lines{ii+kk});
        end
        
    end
    % poles
    if ch1=='P'
        np = str2double(strtrim(ltemp(6:end)));
        p = zeros(np,2);
        for kk=1:np
           %ltemp0 = lines{ii+kk};
           %p(kk,1) = str2double(strtrim(ltemp0(1:14)));
           %p(kk,2) = str2double(strtrim(ltemp0(15:end)));
           p(kk,:) = str2num(lines{ii+kk});
        end
    end
    % constant: CONSTANT = A0 * SENSITIVITY (c = A0*k)
    if ch1=='C'
        c = str2double(strtrim(ltemp(9:end)));
    end
end

% convert to complex
z = complex(z(:,1),z(:,2));
p = complex(p(:,1),p(:,2));

% convert to velocity in m/s by addding a pole (or removing a zero)
if ideriv==0
    disp('read_pzfile.m: default displacement response (m to count)');
elseif ideriv==1
    disp('read_pzfile.m: velocity response (m/s to count)');
    p = [p ; 0]; 
elseif ideriv==2
    disp('read_pzfile.m: acceleration response (m/s^2 to count)');
    p = [p ; 0 ; 0];
end
    
%----------------------
% EXAMPLE

if 0==1
   idir = '/home/admin/databases/SUMATRA/data/CAN_test/';
   sfile = 'SAC_PZs_G_CAN_LHZ__1989.153.00.00.00.0000_2006.344.02.60.60.99999';
   pzfile = [idir sfile];
   [p,z,c,A0,k] = read_pzfile(pzfile,1);
   c, A0*k
end

%==========================================================================

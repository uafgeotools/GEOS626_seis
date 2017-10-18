function y2 = bandpass(t,y,fa,fb)
%BANDPASS apply bandpass filter [fa,fb] on a time series
%
% INPUT
%   t       n x 1 vector of time points (s)
%   y       n x 1 vector of amplitude variations
%   fa      lower frequency limit of bandpass (Hz)
%   fb      upper frequency limit of bandpass (Hz)
%
% OUTPUT
%   y2      n x 1 vector of amplitude variations, bandpass-filtered between
%              frequencies of fa and fb
%
% Code from Muller and MacDonald, Ice Ages and Astronomical Causes, p. 294
%

f1 = min([fa fb]);
f2 = max([fa fb]);

% demean and then take fft
y = y - mean(y);
Hf = fft(y);

% define frequency vector
n = length(t);
dt = t(2) - t(1);
fNyq = 1/(2*dt);
fvec = linspace(0, 2*fNyq, n)';

% find indices that delimit fa,fb,-fa,-fb
[~, k1] = min(abs(fvec - f1));
[~, k2] = min(abs(fvec - f2));
k3 = n - k2 + 2;
k4 = n - k1 + 2;

% zero out entries outside the frequency range [fa,fb]
Hf(1:k1)  = zeros(k1, 1);
Hf(k2:k3) = zeros(k3-k2+1, 1);
Hf(k4:n)  = zeros(k1-1, 1);

% real part of inverse fourier transform
y2 = real(ifft(Hf));

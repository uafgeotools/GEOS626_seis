function [Astack,f,A] = w2fstack(w,f1,f2,n)
%W2STACK compute sectral stack from a vector of waveform objects
%
% This assumes that you have previously computed the spectra using
% wf_fft.compute, which applies the headers fft_freq and fft_amp.
%
% Interpolation is required, since each waveform object may have a
% different frequency vector (for example, if dt differs, then so will df).
%
% INPUT
%   w      vector of waveform objects with precomputed spectra
%   f1,f2  frequency range for stacking, Hz
%   n      number of points to use in frequency domain
%
% OUTPUT
%   Astack n x 1 sum of all spectra
%   f      n x 1 vector of frequencies
%   A      n x nw set of all spectra
%

% compute stack of all spectra
nw = length(w);
f = linspace(f1,f2,n);
A = zeros(n,nw);
for ii=1:nw
    f0 = get(w(ii),'fft_freq');  % Hz
    A0 = get(w(ii),'fft_amp');
    A(:,ii) = interp1(f0,A0,f);
end
Astack = sum(A,2);

function [CSM,freqs,n] = csm(signal,window,noverlap,nfft,fs)

% CSM Compute cross-spectral matrix
%
%function [CSM,freqs] = csm(signal,window,noverlap,nfft,fs)
%
% Input:
%   signal: 2D array, size [N,M] with time data of N samples from M microphones
%   window: Window, e.g. Hann, Hamming, etc. with length = nfft
%   noverlap : Overlap in samples of time segments
%   nfft: Number of samples of each time-segment
%   fs: Sampling frequency
%
% Output: 
%   CSM: Cross-spectral matrix
%   freqs: Frequency vector
%
% Author: Oliver Lylloff
% Date: 8/12/13
% Latest revision: 8/12/13
%
% http://github.com/1oly/mscthesis
% 
% Inspired by FFTANALYSER TOOLBOX by Quentin Leclère
%

[N,M] = size(signal);
n = fix((N-noverlap)/(length(window)-noverlap));    % Number of frames

if round(nfft/2) == nfft/2
    n_spec = nfft/2 + 1;
else
    n_spec = (nfft+1)/2;
end

df = round(fs/nfft);  % Spectral resolution
freqs = (0:n_spec-1)*df;    % Frequency vector

% Initialize
frame = zeros(nfft,M);
spectrum = zeros(n_spec,M);
CSM = zeros(M,M,n_spec);
comp = sum(abs(window).^2)/nfft;
window = window/sqrt(comp);     % Power correction, such that window power = 1
pos = 1;

% First frame
for i = 1:M         
    frame(:,i) = signal(pos:pos+nfft-1,i).*window;   % Current frame of x
end

h = waitbar(0,'Computing CSM');
for i = 1:n
    waitbar(i/n)
    
    % Compute spectrum:
    s = fft(frame)./nfft;
    spectrum(1,:) = s(1,:);    % DC component
    switch fix(nfft/2) == nfft/2
        case 1
            spectrum(2:nfft/2,:) = 2*abs(s(2:nfft/2,:)).*exp(1j*angle(s(2:nfft/2,:)));
            spectrum(nfft/2+1,:) = s(nfft/2+1,:);
        case 0
            spectrum(2:(nfft+1)/2,:) = 2*abs(s(2:(nfft+1)/2,:)).*exp(1j*angle(s(2:(nfft+1)/2,:)));
    end
    
    pos = pos + (nfft - noverlap);  % Update position
    
    % Assemble CSM:
    for j = 1:M
        % Update diagonal of cross-spectral matrix
        CSM(j,j,:) = CSM(j,j,:) + permute(abs(spectrum(:,j)).^2,[2 3 1]);
        % Update lower triangular Rmn
        for k = j+1:M
            temp = CSM(j,k,:) + permute(conj(spectrum(:,j)).*spectrum(:,k),[2 3 1]);
            CSM(j,k,:) = temp;
        end
        % Update new frame
        if i<n
            frame(:,j) = signal(pos:pos+nfft-1,j).*window;
        end
    end
end
close(h);

% Average CSM
CSM = CSM/n/2;      % rms value

% Assemble upper triangular matrix
for j = 1:M
   for k = j+1:M
       CSM(k,j,:) = conj(CSM(j,k,:));
   end
end

end
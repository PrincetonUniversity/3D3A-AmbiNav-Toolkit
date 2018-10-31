function Y = AmbiNav_ForwardSTFT(x, window, noverlap, nfft)

winLen = length(window);
if nargin < 4
    nfft = winLen;
end

Y = spectrogram(x, window, noverlap, nfft, 'twosided'); % NFFT x numTimeFrames

end
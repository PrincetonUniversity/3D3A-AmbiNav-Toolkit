function x = AmbiNav_InverseSTFT(Y, window, noverlap, nfft)

winLen = length(window);
if nargin < 4
    nfft = winLen;
end

invWindow = ones(size(window));
invWindow(window~=0) = 1./window(window~=0);

nCols = size(Y,2);
xLen = (nCols - 1)*(winLen - noverlap) + winLen;
xMat = ifft(Y,nfft,1,'symmetric');

x = zeros(xLen,1);
numFrames = zeros(xLen,1);
for ii = 1:nCols
    indx = (1:winLen) + ((ii-1)*(winLen-noverlap));
    x(indx) = x(indx) + invWindow.*xMat(1:winLen,ii);
    numFrames(indx) = numFrames(indx) + (1 - double(window==0));
end
x(numFrames~=0) = x(numFrames~=0)./numFrames(numFrames~=0);

end
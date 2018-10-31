function ao = thiergart2013(ai, u, Lo, r, Fs, nfft, refMic, force2D)

if nargin < 6
    nfft = [];
end

if nargin < 7
    refMic = 1;
end

if nargin < 8
    force2D = false;
end

% Analyze sound field and construct a model
sf = thiergart2013_analysis(ai, u, Fs, nfft, force2D);

% Synthesize modeled sound field
ao = thiergart2013_synthesis(sf, u, Lo, r, refMic);

end
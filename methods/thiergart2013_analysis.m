function sf = thiergart2013_analysis(ai, u, Fs, nfft, force2D)
%THIERGART2013_ANALYSIS Time-frequency domain analysis of a sound field.
%   SF = THIERGART2013_ANALYSIS(B,U,FS) returns a struct SF containing
%   sound field properties computed via time-frequency domain analysis and
%   source triangulation, given ambisonics signals B measured from
%   positions U and at a sampling rate FS.
%
%   SF = THIERGART2013_ANALYSIS(B,U,FS,NFFT) additionally specifies the
%   NFFT to be used in taking the short-time Fourier transform.
%
%   SF = THIERGART2013_ANALYSIS(B,U,FS,NFFT,FORCE2D) optionally forces all
%   sound field analysis and modeling to be taken in the horizontal plane
%   (i.e., 2D only) by ignoring the third element of each position vector.
%
%   See also THIERGART2013, AMBINAV_TRIANGULATESOURCE.

%   ==============================================================================
%   This file is part of the 3D3A AmbiNav Toolkit.
%   
%   Joseph G. Tylka <josephgt@princeton.edu>
%   3D Audio and Applied Acoustics (3D3A) Laboratory
%   Princeton University, Princeton, New Jersey 08544, USA
%   
%   MIT License
%   
%   Copyright (c) 2018 Princeton University
%   
%   Permission is hereby granted, free of charge, to any person obtaining a copy
%   of this software and associated documentation files (the "Software"), to deal
%   in the Software without restriction, including without limitation the rights
%   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
%   copies of the Software, and to permit persons to whom the Software is
%   furnished to do so, subject to the following conditions:
%   
%   The above copyright notice and this permission notice shall be included in all
%   copies or substantial portions of the Software.
%   
%   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
%   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
%   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
%   SOFTWARE.
%   ==============================================================================

%   References:
%     [1] Thiergart et al. (2013) Geometry-Based Spatial Sound Acquisition
%         Using Distributed Microphone Arrays.

delta = AmbiNav_ArraySpacing(u);
IRLen = size(ai{1},1);

if nargin < 5 || isempty(force2D)
    force2D = false;
end

%% Set STFT parameters
ovlp = 0.5; % assume 50% overlap
if nargin < 4 || isempty(nfft)
    sf.nfft = 2^nextpow2((1/(1-ovlp)) * Fs * delta / getSoundSpeed()); % Eq. (16)
else
    sf.nfft = nfft;
    maxDelta = (getSoundSpeed() / Fs) * sf.nfft * (1-ovlp); % Eq. (16)
    if delta > maxDelta
        warning('Array spacing exceeds maximum for this FFT length; try increasing NFFT.');
    end
end
sf.noverlap = round(ovlp * sf.nfft);
sf.window = hamming(sf.nfft);
if sf.nfft > IRLen
    warning('Due to a large array spacing and short audio length, the STFT has very few time frames.')
end
sf.padSTFT = true; % Zero-pad signal before computing STFT; see getForwardSTFT
sf.numTimeFrames = STFT_len2part(IRLen, sf.nfft, sf.noverlap, sf.padSTFT);
sf.kVec = f2k(getFreqVec(Fs,sf.nfft));
sf.specLen = 1 + sf.nfft/2;

%% Estimate Source DOAs & Reference Pressure
numMics = numel(ai);
if numMics < 2
    error('Not enough microphones!');
end
sf.p = cell(size(ai));
A0_stft = cell(size(ai));
A1_stft = cell(size(ai));
A2_stft = cell(size(ai));
A3_stft = cell(size(ai));
r_I = cell(numMics,sf.numTimeFrames);
s_p = cell(numMics,sf.numTimeFrames,sf.specLen);
sf.psi = cell(numMics,sf.numTimeFrames);
for pp = 1:numMics
    A0_stft{pp} = getForwardSTFT(ai{pp}(:,1), sf.window, sf.noverlap, sf.nfft, sf.padSTFT);
    A1_stft{pp} = getForwardSTFT(ai{pp}(:,2), sf.window, sf.noverlap, sf.nfft, sf.padSTFT);
    A2_stft{pp} = getForwardSTFT(ai{pp}(:,3), sf.window, sf.noverlap, sf.nfft, sf.padSTFT);
    A3_stft{pp} = getForwardSTFT(ai{pp}(:,4), sf.window, sf.noverlap, sf.nfft, sf.padSTFT);
    for ii = 1:sf.numTimeFrames
        A_temp = [A0_stft{pp}(:,ii) A1_stft{pp}(:,ii) A2_stft{pp}(:,ii) A3_stft{pp}(:,ii)];
        [r_I{pp,ii}, sf.psi{pp,ii}] = merimaa2005(A_temp, 'normalize');
        for kk = 1:sf.specLen
            % DOA for mic p at each time-frequency bin
            s_p{pp,ii,kk} = r_I{pp,ii}(kk,:);
        end
    end
    sf.p{pp} = A0_stft{pp} * sqrt(4*pi / ambNormSquared(0,'N3D')); % Pressure at Each Mic
end

%% Estimate Source Positions
sf.s_0 = cell(sf.numTimeFrames,sf.specLen);
for ii = 1:sf.numTimeFrames
    for kk = 1:sf.specLen
        % Triangulated source position for each time-frequency bin
        sf.s_0{ii,kk} = AmbiNav_TriangulateSource(u, s_p(:,ii,kk), force2D);
    end
end

end
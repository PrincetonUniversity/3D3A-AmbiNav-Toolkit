function a_out = thiergart2013_synthesis(sf, u, Lo, r, refMic)
%THIERGART2013_SYNTHESIS Ambisonics rendering of a modeled sound field.
%   A = THIERGART2013_SYNTHESIS(SF,U,L,R) returns the rendered ambisonics
%   signals A, up to order L, at the desired listening position R given a
%   sound field model SF and microphone positions U.
%
%   A = THIERGART2013_SYNTHESIS(SF,U,L,R,REFMIC) additionally specifies
%   which microphone to be used as a reference. By default, the nearest
%   microphone to the listening position is used.
%
%   See also THIERGART2013, THIERGART2013_ANALYSIS.

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

if nargin < 5 || isempty(refMic)
    % By default, use the nearest microphone as a reference
    distVec = zeros(numel(u),1);
    for ii = 1:numel(u)
        distVec(ii) = norm(u{ii} - r);
    end
    refMic = find(distVec == min(distVec),1,'first');
end

%% Estimate Direct and Diffuse Pressures at IPLS
S_dir  = zeros(sf.specLen,sf.numTimeFrames);
S_diff = zeros(sf.specLen,sf.numTimeFrames);
for ii = 1:sf.numTimeFrames
    for kk = 1:sf.specLen
        d_source = sf.s_0{ii,kk} - u{refMic};
        S_ref = sf.p{refMic}(kk,ii); % Pressure at Reference Mic
        
        Gamma_ref = (1/sf.psi{refMic,ii}(kk)) - 1; % Eqs. (26) and (27)
        G_dir_ref = sqrt(Gamma_ref / (1 + Gamma_ref)); % Eq. (8)
        S_dir_ref = G_dir_ref * S_ref; % Eq. (7)
        H_dir = conj(ambRadialFilter(0,sf.kVec(kk),norm(d_source))); % = conj((1/norm(d_source)) * exp(1i*sf.kVec(kk)*norm(d_source)));
            % NOTE: complex conjugatation is needed because the STFT uses
            % the standard FFT convention, whereas the acoustic potential
            % is typically computed using the conjugated Fourier convetion.
        S_dir(kk,ii) = S_dir_ref ./ H_dir; % Eq. (6)
        
        G_diff_ref = sqrt(1 / (1 + Gamma_ref)); % Eq. (14a)
        S_diff_ref = G_diff_ref * S_ref; % Eq. (13)
        S_diff(kk,ii) = S_diff_ref ./ 1; % Eq. (12), with H_diff = 1
    end
end

%% Assemble Output Signals
outLen = STFT_part2len(sf.numTimeFrames, sf.nfft, sf.noverlap, sf.padSTFT);
No = (Lo + 1)^2;
a_out = zeros(outLen,No);
for nn = 1:No
    [l,m] = getAmbOrder(nn-1);
    Q_factor = sqrt(ambNormSquared(l,'N3D') / (4*pi));
    B_n = zeros(sf.nfft,sf.numTimeFrames);
    for ii = 1:sf.numTimeFrames
        for kk = 1:sf.specLen
            d_source = sf.s_0{ii,kk} - r;
            C_factor = ambSphericalHarmonicY(l, m, d_source, 'N3D');
            H_dir = conj(ambRadialFilter(l,sf.kVec(kk),norm(d_source)));
                % See note about complex conjugation above
            B_n(kk,ii) = (C_factor * H_dir * S_dir(kk,ii)) + (Q_factor * 1 * S_diff(kk,ii)); % Eq. (33), with H_diff = 1
        end
    end
    a_out(:,nn) = getInverseSTFT(B_n, sf.window, sf.noverlap, sf.nfft, sf.padSTFT);
end

end
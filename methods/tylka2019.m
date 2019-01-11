function Ao = tylka2019(Ai, u, Lo, r, kVec)
%TYLKA2019 Ambisonics navigation using hybrid interpolation filters.
%   A = TYLKA2019(B,U,L,R,K) computes the interpolated ambisonics signals
%   A, up to order L and interpolated to position vector R (given in
%   Cartesian coordinates), given the ambisonics signals B measured from
%   positions U, and for angular wavenumber K.
%
%   B and U should both be cell arrays with the same number of elements.
%
%   K may be a vector, in which case SIZE(B{1},1) must be LENGTH(K) and A
%   will be LENGTH(K)-by-(LO+1)^2.
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also TYLKA2016, SOUTHERN2009.

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
%     [1] Tylka and Choueiri (2019) A Parametric Method for Virtual
%         Navigation Within an Array of Ambisonics Microphones (Under
%         Review).

Li = sqrt(size(Ai{1},2)) - 1;
xoFreqk = xoFreqModel(u, r, Li);
xoIndx = find(kVec >= xoFreqk,1,'first') - 1;

% Split frequency ranges
AiL = cell(size(Ai));
AiH = cell(size(Ai));
for ii = 1:numel(Ai)
    AiL{ii} = Ai{ii}(1:xoIndx,:);
    AiH{ii} = Ai{ii}((xoIndx+1):end,:);
end

% Low-frequency calculation
AoL = tylka2016(AiL, u, Lo, r, kVec(1:xoIndx));

% High-frequency calcuation
AoH = southern2009(AiH, u, Lo, r);

% Recombine frequency ranges
Ao = cat(1,AoL,AoH);

end

function xoFreqk = xoFreqModel(u, r, Li)
navDist = zeros(size(u));
for ii = 1:numel(u)
    navDist(ii) = norm(u{ii}-r);
end
switch numel(u)
    case 1
        xoFreqk = Li / min(navDist);
%         xoFreqk = 1 / min(navDist);
    case 2
        xoFreqk = Li / AmbiNav_ArraySpacing(u);
%         xoFreqk = 1 / ((min(navDist)*max(navDist)) / AmbiNav_ArraySpacing(u));
    otherwise
        warning('Hybrid crossover frequency is not well-established for P > 2 microphones.');
        xoFreqk = Li / AmbiNav_ArraySpacing(u);
%         xoFreqk = 1 / max(navDist);
end %% TODO: generalizing to P > 2 needs to be investigated more
end
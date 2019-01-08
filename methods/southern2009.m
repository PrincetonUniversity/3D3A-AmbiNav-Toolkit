function Ao = southern2009(Ai, u, Lo, r)
%SOUTHERN2009 Ambisonics navigation using linear interpolation.
%   A = SOUTHERN2009(B,U,L,R) computes the interpolated ambisonics signals
%   A, up to order L and interpolated to position vector R (given in
%   Cartesian coordinates), given the ambisonics signals B measured from
%   positions U.
%
%   B and U should both be cell arrays with the same number of elements.
%
%   See also AMBINAV_INTERPOLATIONWEIGHTS.

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
%     [1] Southern et al. (2009) Rendering walk-through auralisations using
%         wave-based acoustical models.

narginchk(4,4);

if numel(r) ~= 3
    error('Translation vector D should have three elements.');
end

numMics = numel(Ai);
if numel(u) ~= numMics
    error('Microphone signals and positions should have the same number of elements.');
end

No = (Lo + 1)^2;
w = AmbiNav_InterpolationWeights(u,r);

Ao = w(1)*Ai{1}(:,1:No);
for nn = 2:numMics
    Ao = Ao + w(nn)*Ai{nn}(:,1:No);
end

end
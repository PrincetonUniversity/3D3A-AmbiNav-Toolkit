function T = AmbiNav_PlaneWaveTranslation(rp, d, kVec)
%AMBINAV_PLANEWAVETRANSLATION Plane-wave translation coefficients matrix.
%   T = AMBINAV_PLANEWAVETRANSLATION(RP,D,K) computes the plane-wave
%   translation coefficients matrix T, for a grid of plane-wave directions
%   RP (given in Cartesian coordinates), translation position vector D, and
%   for angular wavenumber K.
%
%   K may be a vector, in which case T is LENGTH(K)-by-SIZE(RP,1).
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also SCHULTZ2013.

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

narginchk(3,3);

if numel(d) ~= 3
    error('Translation vector D should have three elements.');
end

rpHat = normalizeVector(rp, 2);
T = exp(1i*shiftdim(kVec) * (rpHat * d.').');

end
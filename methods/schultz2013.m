function [Ao, MUo] = schultz2013(Ai, Lo, pwGrid, d, kVec, wQList)
%SCHULTZ2009 Ambisonics navigation using plane-wave translation.
%   B = SCHULTZ2013(A,LO,RP,D,K) computes the translated ambisonics
%   potentials B, up to order LO, given the ambisonics potentials A, a grid
%   of plane-wave directions RP, a translation position vector D (given in
%   Cartesian coordinates), and for angular wavenumber K. K may be a
%   vector, in which case SIZE(A,1) must be LENGTH(K) and B is
%   LENGTH(K)-by-(LO+1)^2. The N3D ambisonics normalization convention is
%   assumed.
%
%   B = SCHULTZ2013(A,LO,RP,D,K,WQ) uses quadrature weights WQ. If
%   unspecified, WQ is given by 1/SIZE(RP,1).
%
%   See also AMBINAV_PLANEWAVETRANSLATION.

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
%     [1] Schultz and Spors (2013)

narginchk(5,6);

if nargin < 6 || isempty(wQList)
    wQList = 1/size(pwGrid,1);
end

if numel(d) ~= 3
    error('Translation vector D should have three elements.');
end

Ni = size(Ai,2);

MUi = Ai*AmbiNav_sphericalHarmonicY(sqrt(Ni)-1, pwGrid);
T = AmbiNav_planeWaveTranslation(-pwGrid, d, kVec);
MUo = (MUi.*T)*diag(wQList);
Ao = MUo*(AmbiNav_sphericalHarmonicY(Lo, pwGrid).');

end
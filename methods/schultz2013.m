function [Ao, MUo] = schultz2013(Ai, Lo, pwGrid, d, kVec, wQList, pinvFlag)
%SCHULTZ2009 Ambisonics navigation using plane-wave translation.
%   B = SCHULTZ2013(A,LO,RP,D,K) computes the translated ambisonics
%   potentials B, up to order LO, given the ambisonics potentials A, a grid
%   of plane-wave source directions RP (given in Cartesian coordinates), a
%   translation position vector D, and for angular wavenumber K.
%
%   K may be a vector, in which case SIZE(A,1) must be LENGTH(K) and B will
%   be LENGTH(K)-by-(LO+1)^2.
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   B = SCHULTZ2013(A,LO,RP,D,K,WQ) uses quadrature weights WQ. If
%   unspecified, WQ is given by 1/SIZE(RP,1).
%
%   B = SCHULTZ2013(A,LO,RP,D,K,WQ,PINVFLAG) additionally specifies the 
%   plane-wave conversion method to use, where if PINVFLAG evaluates to
%   true, a least-squares inversion is taken to compute the plane-wave
%   signals. By default, PINVFLAG is false.
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
%     [1] Schultz and Spors (2013) Data-based Binaural Synthesis Including 
%         Rotational and Translatory Head-Movements.

narginchk(5,7);

if numel(d) ~= 3
    error('Translation vector D should have three elements.');
end

if nargin < 6 || isempty(wQList)
    wQList = [];
end

if nargin < 7 || isempty(pinvFlag)
    pinvFlag = false;
end

MUi = a2mu(Ai,pwGrid,'N3D',pinvFlag,wQList);
T = AmbiNav_PlaneWaveTranslation(pwGrid, d, kVec);
MUo = MUi.*T;
Ao = mu2a(MUo,pwGrid,Lo,wQList,'N3D');

end
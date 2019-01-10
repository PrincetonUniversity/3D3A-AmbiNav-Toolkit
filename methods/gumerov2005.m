function Ao = gumerov2005(Ai, Lo, d, kVec)
%GUMEROV2005 Ambisonics navigation using translation coefficients.
%   A = GUMEROV2005(B,L,D,K) computes the translated ambisonics potentials
%   A, up to order L, given the ambisonics potentials B, a translation
%   position vector D (given in Cartesian coordinates), and for angular
%   wavenumber K.
%
%   K may be a vector, in which case SIZE(B,1) must be LENGTH(K) and A will
%   be LENGTH(K)-by-(LO+1)^2.
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also AMBINAV_TRANSLATION.

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
%     [1] Gumerov and Duraiswami (2005) Fast Multipole Methods for the
%         Helmholtz Equation in Three Dimensions.
%     [2] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

narginchk(4,4);

if numel(d) ~= 3
    error('Translation vector D should have three elements.');
end

kLen = length(kVec);

Ni = size(Ai,2);
No = (Lo + 1)^2;

T = AmbiNav_Translation(sqrt(Ni)-1, Lo, d, kVec);

Ao = zeros(kLen,No);
% for kk = 1:kLen
%     Ao(kk,:) = Ai(kk,:) * T(:,:,kk);
% end
for nn = 1:No
    Ao(:,nn) = sum(Ai .* squeeze(T(:,nn,:)).',2);
end

end
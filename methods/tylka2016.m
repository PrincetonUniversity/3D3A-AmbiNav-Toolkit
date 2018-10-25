function Ao = tylka2016(Ai, u, Lo, d, kVec, beta)
%TYLKA2016 Ambisonics navigation using least-squares interpolation filters.
%   B = TYLKA2016(A,U,LO,D,K) computes the interpolated ambisonics signals
%   B, up to order LO and interpolated to position vector D (given in
%   Cartesian coordinates), given the ambisonics signals A measured from
%   positions U, and for angular wavenumber K.
%
%   A and U should both be cell arrays with the same number of elements.
%
%   K may be a vector, in which case SIZE(A{1},1) must be LENGTH(K) and B
%   will be LENGTH(K)-by-(LO+1)^2.
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also AMBINAV_INTERPOLATIONFILTERS.

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
%     [1] Tylka and Choueiri (2016) Soundfield Navigation using an Array of
%         Higher-Order Ambisonics Microphones.

narginchk(5,6);

if nargin < 6
    beta = []; % Use default regularization function
end

if numel(d) ~= 3
    error('Translation vector D should have three elements.');
end

numMics = numel(Ai);
if numel(u) ~= numMics
    error('Microphone signals and positions should have the same number of elements.');
end

Li = sqrt(size(Ai{1},2)) - 1;
kLen = length(kVec);

Acat = cell2mat(Ai(:).');
Ao = zeros(kLen,(Lo + 1)^2);
F = AmbiNav_InterpolationFilters(Li,Lo,u,d,kVec,beta);
for kk = 1:length(kVec)
    Ao(kk,:) = Acat(kk,:) * F(:,:,kk).';
end

end
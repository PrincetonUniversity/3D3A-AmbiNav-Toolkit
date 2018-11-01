function Y = AmbiNav_SphericalHarmonicY(L, R)
%AMBINAV_SPHERICALHARMONICY Real-valued spherical harmonic function.
%   Y = AMBINAV_SPHERICALHARMONICY(L,R) computes the real-valued, N3D
%       normalized spherical harmonics, up to order L and for positions R.
%       The N3D spherical harmonic convention is described by Zotter [1].
%
%   Note:
%       L must be a scalar.
%       R may be a P-by-3 matrix of directions, where each row is a
%           Cartesian vector.
%       Y will be a (L + 1)^2-by-P matrix.

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
%     [1] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

% Needs at least 2 input arguments
if nargin < 2
    error('Not enough input arguments.');
end

% Compute spherical harmonic matrix
Y = zeros((L + 1)^2, size(R,1));
for l = 0:L
    for m = -l:l
        acn = l*(l + 1) + m;
        Y(acn + 1,:) = ambSphericalHarmonicY(l, m, R, 'N3D');
    end
end

end
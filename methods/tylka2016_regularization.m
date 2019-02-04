function beta = tylka2016_regularization(k, u, r, L, TYPE)
%TYLKA2016_REGULARIZATION Regularization function for pinv interpolation.
%   BETA = TYLKA2016_REGULARIZATION(K,U,R) computes a regularization
%   function BETA for computing regularized least-squares interpolation
%   filters from microphones at positions U to a listener at position R,
%   for angular wavenumber K.
%
%   K may be a vector, in which case B is LENGTH(K)-by-1.
%
%   BETA = TYLKA2016_REGULARIZATION(K,U,R,L) additionally specifies the
%   input ambisonics order L (by default, L = 1), which may be used by some
%   regularization functions.
%
%   BETA = TYLKA2016_REGULARIZATION(K,U,R,L,TYPE) additionally specifies
%   the TYPE of regularization function. Currently implemented options are
%   'tylka2016' (default) and 'tylka2019'.
%
%   See also AMBINAV_INTERPOLATIONFILTERS, TYLKA2016, TYLKA2019.

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
%     [2] Tylka and Choueiri (2019) A Parametric Method for Virtual
%         Navigation Within an Array of Ambisonics Microphones (Under
%         Review).

narginchk(3,5);

if nargin < 4 || isempty(L)
    L = 1;
end

if nargin < 5 || isempty(TYPE)
    TYPE = 'tylka2016';
end

beta = zeros(length(k),1);

numMics = numel(u);
navDist = zeros(numMics,1);
for ll = 1:numMics
    navDist(ll) = norm(r - u{ll});
end

switch lower(TYPE)
    case 'tylka2016'
        k0 = 1 / AmbiNav_ArraySpacing(u);
    case 'tylka2019'
        switch numMics
            case 1
                k0 = 1 / min(navDist);
            case 2
                k0 = AmbiNav_ArraySpacing(u) / (min(navDist) * max(navDist));
            otherwise
                k0 = 1 / max(navDist);
        end
    otherwise
        error('No known regularization type %s.', TYPE);
end

Gpi = db2mag(30);
beta(k~=0) = abs((Gpi * 1i * (k(k~=0) / k0) + 1) ./ (1i * (k(k~=0) / k0) + Gpi));

end
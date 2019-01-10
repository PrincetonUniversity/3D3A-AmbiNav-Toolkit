function F = AmbiNav_InterpolationFilters(Li, Lo, u, r, kVec, beta)
%AMBINAV_INTERPOLATIONFILTERS Regularized least-squares interpolation filters.
%   F = AMBINAV_INTERPOLATIONFILTERS(LI,LO,U,R,K) computes least-squares
%   ambisonics interpolation filters F, up to order LO for interpolation to
%   position vector R (given in Cartesian coordinates), given P measurement
%   positions U.
%
%   K may be a vector, in which case F is (LO+1)^2-by-(P*(LI+1)^2)-by-LENGTH(K).
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also AMBINAV_INTERPOLATIONWEIGHTS, AMBINAV_TRANSLATION, TYLKA2016.

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

if numel(r) ~= 3
    error('Position vector R should have three elements.');
end

weights = AmbiNav_InterpolationWeights(u,r);
weights(weights < max(weights)/1e10) = 0;

% Default regularization function
if nargin < 6 || isempty(beta)
    beta = regFunction(kVec,u(weights~=0),r);
end

numMics = numel(u);
Ni = (Li + 1)^2;
No = (Lo + 1)^2;
Lmax = floor(sqrt(Ni * numMics)) - 1;
Nmax = (Lmax + 1)^2;
kLen = length(kVec);

indx = cell(numMics,1);
T_cell = cell(numMics,1);
for ll = 1:numMics
    indx{ll} = 1+(ll-1)*Ni:ll*Ni;
    
    % Each cell contains a matrix of size Nmax-by-Ni-by-K
    if weights(ll) ~= 0
        T_cell{ll} = AmbiNav_Translation(Lmax, Li, u{ll} - r, kVec);
    else
        T_cell{ll} = zeros(Nmax, Ni, kLen);
    end
end

F = zeros(No, numMics * Ni, kLen);
W = zeros(numMics * Ni);
M = zeros(Ni * numMics, Nmax);
for kk = 1:kLen
    if kVec(kk) == 0
        for ll = 1:numMics
            F(:,indx{ll},kk) = weights(ll)*eye(No,Ni);
        end
    else
        for ll = 1:numMics
            M(indx{ll},:) = sqrt(weights(ll))*(T_cell{ll}(:,:,kk).');
            W(indx{ll},indx{ll}) = sqrt(weights(ll))*eye(Ni);
        end
        
        while any(any(isnan(pinv(M))))
            weights = ditherWeights(weights);
            for ll = 1:numMics
                M(indx{ll},:) = sqrt(weights(ll))*(T_cell{ll}(:,:,kk).');
            end
        end
        
        [U, S, V] = svd(M);
        TH = zeros(size(S,2));
        sigma = diag(S);
        numSigma = length(sigma);
        
        TH(1:numSigma,1:numSigma) = diag(sigma.^2./(sigma.^2 + ((max(sigma)/1000) * beta(kk))));
        F(:,:,kk) = V(1:No,:)*TH*pinv(S)*U'*W;
    end
end

end

function w = ditherWeights(w)

d = 0.01*min(min(w),1-max(w));
w(w==max(w)) = w(w==max(w)) + d;
w(w==min(w)) = w(w==min(w)) - d;

end

function beta = regFunction(k, u, r)

beta = zeros(length(k),1);

numMics = numel(u);
navDist = zeros(numMics,1);
for ll = 1:numMics
    navDist(ll) = norm(r - u{ll});
end

switch numMics
    case 1
        k0 = 1 / navDist;
    case 2
        k0 = norm(u{1} - u{2}) / (min(navDist) * max(navDist));
    otherwise
        delta = AmbiNav_ArraySpacing(u);
        k0 = 1 / delta;
end

Gpi = db2mag(30);
beta(k~=0) = abs((Gpi * 1i * (k(k~=0) / k0) + 1) ./ (1i * (k(k~=0) / k0) + Gpi));

end
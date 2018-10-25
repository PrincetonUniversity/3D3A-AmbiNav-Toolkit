function W = AmbiNav_InterpolationWeights(X, Xq)
%AMBINAV_INTERPOLATIONWEIGHTS Linear interpolation weights.
%   W = AMBINAV_INTERPOLATIONWEIGHTS(X,XQ) computes linear interpolation
%   weights W for a function F measured at positions X and for query points
%   XQ. X should be either an N-by-D matrix or a cell array with N
%   elements, each of which should be a 1-by-D vector, where, N is the
%   number of measurement positions and D is the dimension of the
%   coordinate system. XQ should be either a Q-by-D matrix or a cell array
%   with Q elements, each of which should be a 1-by-D vector, where Q is
%   the number query points. W will then be a Q-by-N matrix.

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

if iscell(X)
    Y = cell2mat(X(:));
else
    Y = X;
end

if iscell(Xq)
    Yq = cell2mat(Xq(:));
else
    Yq = Xq;
end

[N, D] = size(Y);
if size(Yq,2) ~= D
    error('Vectors must be the same length.');
end

% Remove constant dimensions
dimMask = true(1,D);
for dd = 1:D
    dimMask(dd) = ~all(Y(:,dd)==Y(1,dd));
end

Q = size(Yq,1);

W = cat(2,Yq(:,dimMask),ones(Q,1)) / cat(2,Y(:,dimMask),ones(N,1));
% NOTE: this approach can have problems, especially when X has collinear
% points. Need to look into using rank calculation and nearest neighbors to
% pick out only the relevant grid points

% W = zeros(Q,N);
% for ii = 1:N
%     v = zeros(N,1);
%     v(ii) = 1;
%     F = scatteredInterpolant(Y,v,'linear');
%     W(:,ii) = F(Yq);
% end

% Remove negative weights for each query point
if any(W(:)<0)
    for qq = 1:Q
        posMask = W(qq,:)>=0;
        if ~all(posMask)
            W(qq,~posMask) = 0;
            W(qq,posMask) = cat(2,Yq(qq,dimMask),ones(1,1)) / cat(2,Y(posMask,dimMask),ones(sum(posMask),1));
        end
    end
end

end
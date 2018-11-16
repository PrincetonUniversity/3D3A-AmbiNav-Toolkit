function [B, v] = AmbiNav_CheckValidity(A, u, r, s)
%AMBINAV_CHECKVALIDITY Identify valid microphones.
%   [B,V] = AMBINAV_CHECKVALIDITY(A,U,R,S) given a set of ambisonics
%   microphones with signals A and positions U, this function returns the
%   signals B and positions V for a subset of those microphones, each of 
%   which provides a valid description of the sound field at the listening
%   position R given source position(s) S.
%
%   If S is empty or omitted, the position(s) of any sources will be
%   estimated via triangulation.

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

if nargin < 4 || isempty(s)
    d = cell(size(u));
    for ii = 1:numel(u)
        d{ii} = estimateSourceDirection(A{ii}, 'N3D');
    end
    s = {AmbiNav_TriangulateSource(u, d)};
end

if ~iscell(s)
    s = {s};
end

navDist = zeros(size(u));
for ii = 1:numel(u)
    navDist(ii) = norm(u{ii}-r);
end

mask = true(size(u));

for nn = 1:numel(s)
    srcDist = zeros(size(u));
    for ii = 1:numel(u)
        srcDist(ii) = norm(u{ii}-s{nn});
    end
    
    mask = mask & (navDist < srcDist);
end

if ~any(mask(:))
    warning('No valid microphones, using nearest.');
    v = cell(1);
    [v{1}, indx] = findNearest(cell2mat(u(:)), r, 'l2', 1);
else
    indx = find(mask);
    v = u(indx);
end

B = A(indx);

end
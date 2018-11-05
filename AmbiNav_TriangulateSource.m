function s0 = AmbiNav_TriangulateSource(u, s, force2D)
%AMBINAV_TRIANGULATESOURCE Estimated source position via triangulation.
%   S0 = AMBINAV_TRIANGULATESOURCE(U,S) esimates the source position S0
%   given a cell array of microphone positions U and a cell array of
%   corresponding directions-of-arrival S.
%
%   S0 = AMBINAV_TRIANGULATESOURCE(U,S,FORCE2D) optionally forces
%   triangulation to be done in 2D by ignoring the third element of each
%   vector.

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

narginchk(2,3);

s0 = zeros(1,3);

if force2D
    uMat = cell2mat(u(:));
    sMat = cell2mat(s(:));
    u = mat2cell(uMat(:,1:2),ones(size(uMat,1),1));
    s = mat2cell(sMat(:,1:2),ones(size(sMat,1),1));
end

switch numel(u)
    case 2
        if force2D
            [sout, ~] = triangulate_2d(u{1}, u{2}, s{1}, s{2});
        else
            [sout, ~] = triangulate_3d(u{1}, u{2}, s{1}, s{2});
        end
    otherwise
        [sout, ~] = triangulate_opt(u, s);
end

if force2D
    s0(1:2) = sout;
else
    s0 = sout;
end

end

function [s0, c0] = triangulate_2d(u1, u2, s1, s2)

s0 = zeros(size(s1));

x1 = u1(1);
y1 = u1(2);
x2 = u2(1);
y2 = u2(2);
a1 = s1(1);
b1 = s1(2);
a2 = s2(1);
b2 = s2(2);

c0 = [a1, -a2; b1, -b2]\[x2-x1; y2-y1];
s0(1:2) = u1(1:2) + c0(1)*s1(1:2);

end

function [s0, c0] = triangulate_3d(u1, u2, s1, s2)

x1 = u1(1);
y1 = u1(2);
z1 = u1(3);
x2 = u2(1);
y2 = u2(2);
z2 = u2(3);
a1 = s1(1);
b1 = s1(2);
c1 = s1(3);
a2 = s2(1);
b2 = s2(2);
c2 = s2(3);

c0 = [a1, -a2; b1, -b2; c1, -c2]\[x2-x1; y2-y1; z2-z1];
s0 = ((u1 + c0(1)*s1) + (u2 + c0(2)*s2))/2;

end

function [s0, c0] = triangulate_opt(u, s)

nPts = numel(u);
if nPts < 2
    error('Need at least 2 measurement points!');
end

nDims = length(u{1});
switch nDims
    case 2
        [si, ~] = triangulate_2d(u{1}, u{2}, s{1}, s{2});
    case 3
        [si, ~] = triangulate_3d(u{1}, u{2}, s{1}, s{2});
end
ci = triangulation_coeffs(u, s, si);

c0 = fminunc(@(c) triangulation_error(u, s, c), ci);
s0 = average_direction(u, s, c0);

end

function cost = triangulation_error(u, s, c)

x = average_direction(u, s, c);
nPts = numel(u);
costs = zeros(nPts,1);
for mm = 1:nPts
    costs(mm) = norm(x - (u{mm} + c(mm)*s{mm}));
end
cost = sum(costs);

end

function c = triangulation_coeffs(u, s, x)

nPts = numel(u);
c = zeros(nPts,1);
for mm = 1:nPts
    c(mm) = (x-u{mm})/s{mm};
end

end

function x = average_direction(u, s, c)

nPts = numel(u);
x = (u{1} + c(1)*s{1})/nPts;
for mm = 2:nPts
    x = x + (u{mm} + c(mm)*s{mm})/nPts;
end

end
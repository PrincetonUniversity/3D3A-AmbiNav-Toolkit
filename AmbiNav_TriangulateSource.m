function s0 = AmbiNav_TriangulateSource(u, s, force2D)

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
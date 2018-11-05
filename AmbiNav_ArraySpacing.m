function delta = AmbiNav_ArraySpacing(u)
%AMBINAV_ARRAYSPACING Spacing between microphone positions.
%   D = AMBINAV_ARRAYSPACING(U) returns the array spacing D given a cell
%   array of microphone positions U.
%
%   For a single microphone (NUMEL(U) = 1), D = NORM(U{1}).
%
%   For two microphones (NUMEL(U) = 2), D = NORM(U{1} - U{2}).
%
%   Otherwise, D will be the most common element of the set of all pairwise
%   distances between each unique pair of elements of U.

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

numMics = numel(u);
switch numMics
    case 1
        delta = norm(u{1});
    case 2
        delta = norm(u{1} - u{2});
    otherwise
        micDist = zeros(sum(1:(numMics-1)),1);
        nn = 1;
        for ii = 1:numMics
            for jj = (ii+1):numMics
                micDist(nn,1) = norm(u{ii} - u{jj});
                nn = nn + 1;
            end
        end
        delta = mode(micDist); % Approximate array spacing
end

end
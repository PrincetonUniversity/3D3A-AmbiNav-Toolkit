function ao = thiergart2013(ai, u, Lo, r, Fs, nfft, refMic, force2D)
%THIERGART2013 Ambisonics navigation via sound field analysis and modeling.
%   A = THIERGART2013(B,U,L,R,FS) computes the ambisonics signals A, up to
%   order L, at the desired listening position R given ambisonics input
%   signals B recorded at positions U and at sampling rate FS.
%
%   A = THIERGART2013(B,U,L,R,FS,NFFT) additionally specifies the NFFT to
%   be used in taking the short-time Fourier transform.
%
%   A = THIERGART2013(B,U,L,R,FS,NFFT,REFMIC) additionally specifies which
%   microphone to be used as a reference. By default, the nearest
%   microphone to the listening position is used.
%
%   A = THIERGART2013(B,U,L,R,FS,NFFT,REFMIC,FORCE2D) optionally forces
%   all sound field analysis and modeling to be taken in the horizontal
%   plane (i.e., 2D only) by ignoring the third element of each position
%   vector. By default, FORCE2D evaluates to false.
%
%   See also THIERGART2013_ANALYSIS, THIERGART2013_SYNTHESIS.

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
%     [1] Thiergart et al. (2013) Geometry-Based Spatial Sound Acquisition
%         Using Distributed Microphone Arrays.

if nargin < 6
    nfft = [];
end

if nargin < 7
    refMic = [];
end

if nargin < 8
    force2D = false;
end

% Analyze sound field and construct a model
sf = thiergart2013_analysis(ai, u, Fs, nfft, force2D);

% Synthesize modeled sound field
ao = thiergart2013_synthesis(sf, u, Lo, r, refMic);

end
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

AmbiNav_Start;

%% Examples

Fs = 48000;
FFTLen = 2048;
kLen = 1 + FFTLen/2;

freqVec = getFreqVec(Fs,FFTLen);
kVec = f2k(freqVec(1:kLen));
timeVec = getTimeVec(Fs,FFTLen);

a1 = zeros(2048,16); a1(1,1) = 1;
a2 = zeros(2048,16); a2(1,2) = 1;
a3 = zeros(2048,16); a3(1,3) = 1;
a4 = zeros(2048,16); a4(1,4) = 1;
A1 = getPotential(a1,FFTLen,1);
A2 = getPotential(a2,FFTLen,1);
A3 = getPotential(a3,FFTLen,1);
A4 = getPotential(a4,FFTLen,1);
Ain = {A1(1:kLen,:),A2(1:kLen,:),A3(1:kLen,:),A4(1:kLen,:)};

u = {[1 1 0],[1 -1 0],[-1 1 0],[-1 -1 0]};
r = [0.1 0.2 0];

Lout = 1;
Nout = (Lout + 1)^2;

[posMat, wQList] = loadGridFile('fliege_36');

figure()
hold all

%%

Aout = gumerov2005(Ain{1}, Lout, r - u{1}, kVec);
Aout = cat(1,Aout,zeros(FFTLen-kLen,Nout));
aout = getPressure(Aout,FFTLen,1,'symmetric');

plot(timeVec,aout(:,1))

%%

[Aout, ~] = schultz2013(Ain{1}, Lout, posMat, r - u{1}, kVec, wQList);
Aout = cat(1,Aout,zeros(FFTLen-kLen,Nout));
aout = getPressure(Aout,FFTLen,1,'symmetric');

plot(timeVec,aout(:,1))

%%

Aout = southern2009(Ain, u, Lout, r);
Aout = cat(1,Aout,zeros(FFTLen-kLen,Nout));
aout = getPressure(Aout,FFTLen,1,'symmetric');

plot(timeVec,aout(:,1))

%%

Aout = tylka2016(Ain, u, Lout, r, kVec);
Aout = cat(1,Aout,zeros(FFTLen-kLen,Nout));
aout = getPressure(Aout,FFTLen,1,'symmetric');

plot(timeVec,aout(:,1))
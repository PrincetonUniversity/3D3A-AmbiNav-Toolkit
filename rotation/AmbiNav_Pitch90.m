function Qp = AmbiNav_Pitch90(L)
%AMBINAV_PITCH90 Ambisonics rotation of 90 degrees pitch.
%   Q = AMBINAV_PITCH90(L) computes the ambisonic rotation coefficients
%   matrix Q, up to ambisonics order L, for a rotation of 90 degrees pitch.
%
%   See also AMBINAV_COEFFICIENTA, AMBINAV_COEFFICIENTB.

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
%     [1] Gumerov and Duraiswami (2005) Fast Multipole Methods for the
%         Helmholtz Equation in Three Dimensions.
%     [2] Zotter (2009) Analysis and Synthesis of Sound-Radiation with
%         Spherical Arrays.

N = (L + 1)^2;
[Lvec, Mvec] = getAmbOrder(0:N-1);

Q = zeros((2*L + 1)^2);

% Step 1
for ll = 0:(2*L)
    Pl = legendre(ll,0);
    for mi = 0:ll
        % Eq. 189 [2], mo = 0
        Q(getACN(ll,0)+1,getACN(ll,mi)+1) = ((-1)^abs(mi)) * sqrt(2-(~abs(mi))) ...
            * sqrt(factorial(ll-abs(mi))/factorial(ll+abs(mi))) * Pl(abs(mi)+1);
    end
end

% Step 2
for mo = 1:L
    for ll = mo:(2*L - mo)
        coeff = sqrt(2-(~mo)) / (2*AmbiNav_CoefficientB(ll+1,mo-1)*sqrt(2-(~(mo-1))));
        for mi = mo:ll
            % Eq. 190 [2]
            temp1 = (AmbiNav_CoefficientB(ll+1,mi-1)/sqrt(2-(~(mi-1)))) * Q(getACN(ll+1,mo-1)+1,getACN(ll+1,mi-1)+1);
            temp2 = (AmbiNav_CoefficientB(ll+1,-mi-1)/sqrt(2-(~(mi+1)))) * Q(getACN(ll+1,mo-1)+1,getACN(ll+1,mi+1)+1);
            temp = sqrt(2-(~mi))*(temp1 - temp2);
            Q(getACN(ll,mo)+1,getACN(ll,mi)+1) = coeff*(temp ...
                + 2*AmbiNav_CoefficientA(ll,mi)*Q(getACN(ll+1,mo-1)+1,getACN(ll+1,mi)+1));
        end
    end
end

% Step 3
for ll = 1:L
    for mo = 1:ll
        for mi = 0:(mo-1)
            % Eq. 191 [2]
            Q(getACN(ll,mo)+1,getACN(ll,mi)+1) = ((-1)^(mo+mi))*Q(getACN(ll,mi)+1,getACN(ll,mo)+1);
        end
    end
end

% Step 4
for ll = 1:L
    for mo = 1:ll
        for mi = 1:ll
            % Eqs. 176 and 191 [2]
            Q(getACN(ll,-mo)+1,getACN(ll,-mi)+1) = Q(getACN(ll,mo)+1,getACN(ll,mi)+1);
        end
    end
end

% Step 5
Qmask = false(N);
for no = 1:N
    lo = Lvec(no);
    mo = Mvec(no);
    for ni = 1:N
        li = Lvec(ni);
        mi = Mvec(ni);
        if lo == li
            % Eq. 187 [2]
            if (mo >= 0) && (mi >= 0)
                Qmask(no,ni) = ~mod(lo + mo + mi, 2);
            elseif (mo < 0) && (mi < 0)
                Qmask(no,ni) = ~mod(lo + mo + mi + 1, 2);
            end
        end
    end
end

Qp = Qmask .* Q(1:N,1:N);

% TESTED 8 JAN 2019 - PASSED %
% matches Kronlachner's ambix Rotator plug-in up to L = 3 %

end
function Tz = AmbiNav_ZTranslation(kd, L)
%AMBINAV_ZTRANSLATION Ambisonics translation along the z axis.
%   T = AMBINAV_ZTRANSLATION(KD,L) computes the ambisonic translation
%   coefficients matrix T, up to ambisonics order L and for non-dimensional
%   frequency KD, given by product of the angular wavenumber K and the
%   translation distance D.
%
%   KD may be a vector, in which case T is (L+1)^2-by-(L+1)^2-by-LENGTH(KD).
%
%   The ACN/N3D ambisonics normalization convention is assumed.
%
%   See also AMBINAV_TRANSLATION.

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

narginchk(2,2);

kdLen = length(kd);
N = (L + 1)^2;
Tz = zeros(N, (2*L+1)^2, kdLen);
zkd = (kd==0) | (kd < AmbiNav_KDThreshold());
nzkd = ~zkd;

zkdPos = find(zkd);
for iii = 1:sum(zkd)
    Tz(:,:,zkdPos(iii)) = eye(N, (2*L+1)^2);
end

% Step 1
for lp = 0:2*L
    % Eq. 166 [2]; Eq. 3.2.103 [1]
    Tz(getACN(0,0)+1,getACN(lp,0)+1,nzkd) = ((-1)^lp)*sqrt(2*lp+1)*sphericalBesselJ(lp,kd(nzkd));
end

% Step 2
for ll = 1:L
    mm = ll;
    for lp = ll:(2*L-ll)
        % Eq. 163 [2]; Eq. 3.2.104 [1]
        term1 = AmbiNav_CoefficientB(lp+1,mm-1) * Tz(getACN(ll-1,mm-1)+1,getACN(lp+1,mm-1)+1,nzkd);
        term2 = AmbiNav_CoefficientB(lp,-mm) * Tz(getACN(ll-1,mm-1)+1,getACN(lp-1,mm-1)+1,nzkd);
        Tz(getACN(ll,mm)+1,getACN(lp,mm)+1,nzkd) = (-term1 + term2) / AmbiNav_CoefficientB(ll,-mm);
    end
end

% Step 3
for mm = 0:(L-1)
    for ll = (mm+1):L
        for lp = ll:(2*L - ll)
            % Eq. 163 [2]; Eq. 3.2.90 [1]
            term1 = AmbiNav_CoefficientA(lp,mm) * Tz(getACN(ll-1,mm)+1,getACN(lp+1,mm)+1,nzkd);
            term2 = AmbiNav_CoefficientA(lp-1,mm) * Tz(getACN(ll-1,mm)+1,getACN(lp-1,mm)+1,nzkd);
            term3 = AmbiNav_CoefficientA(ll-2,mm);
            if term3 ~= 0
                term3 = term3 * Tz(getACN(ll-2,mm)+1,getACN(lp,mm)+1,nzkd);
            end
            Tz(getACN(ll,mm)+1,getACN(lp,mm)+1,nzkd) = (-term1 + term2 + term3) / AmbiNav_CoefficientA(ll-1,mm);
        end
    end
end

% Step 4
for ll = 1:L
    for lp = ll:L
        for mm = 1:ll
            % Eq. 161 [2]; Eq. 3.2.92 [1]
            Tz(getACN(ll,-mm)+1,getACN(lp,-mm)+1,nzkd) = Tz(getACN(ll,mm)+1,getACN(lp,mm)+1,nzkd);
        end
    end
end

% Step 5
for lp = 0:(L-1)
    for ll = (lp+1):L
        for mm = -lp:lp
            % Eq. 162 [2]; Eq. 3.2.96 [1]
            coeff = (-1)^(ll+lp);
            Tz(getACN(ll,mm)+1,getACN(lp,mm)+1,nzkd) = coeff * Tz(getACN(lp,mm)+1,getACN(ll,mm)+1,nzkd);
        end
    end
end

% Make square matrix
Tz = Tz(:,1:N,:);

% Basis function correction (since we factor out (-i)^l)
[Lvec, ~] = getAmbOrder(0:N-1);
for nn = 1:N
    ll = Lvec(nn);
    for np = 1:N
        lp = Lvec(np);
        if any(Tz(nn,np,nzkd))
            coeff = (-1i)^(ll-lp);
            Tz(nn,np,nzkd) = coeff * Tz(nn,np,nzkd);
        end
    end
end

% Transpose to translate ambisonics signals (rather than basis functions)
for kk = 1:kdLen
    Tz(:,:,kk) = Tz(:,:,kk).';
end

end
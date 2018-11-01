function sf = thiergart2013_analysis(ai, u, Fs, nfft, force2D)

delta = AmbiNav_ArraySpacing(u);
IRLen = size(ai{1},1);

%% Set STFT parameters
ovlp = 0.5; % assume 50% overlap
if nargin < 4 || isempty(nfft)
    sf.nfft = 2^nextpow2((1/(1-ovlp)) * Fs * delta / getSoundSpeed()); % Eq. (16)
else
    sf.nfft = nfft;
    maxDelta = (getSoundSpeed() / Fs) * sf.nfft * (1-ovlp); % Eq. (16)
    if delta > maxDelta
        warning('Array spacing exceeds maximum for this FFT length; try increasing NFFT.');
    end
end
sf.noverlap = round(ovlp * sf.nfft);
sf.window = hamming(sf.nfft);
if sf.nfft > IRLen
    warning('Due to a large array spacing and short audio length, the STFT has very few time frames.')
end
sf.numTimeFrames = fix((IRLen-sf.noverlap)/(sf.nfft-sf.noverlap));
sf.kVec = f2k(getFreqVec(Fs,sf.nfft));
sf.specLen = 1 + sf.nfft/2;

%% Estimate Source DOAs & Reference Pressure
numMics = numel(ai);
if numMics < 2
    error('Not enough microphones!');
end
sf.p = cell(size(ai));
A0_stft = cell(size(ai));
A1_stft = cell(size(ai));
A2_stft = cell(size(ai));
A3_stft = cell(size(ai));
r_I = cell(numMics,sf.numTimeFrames);
s_p = cell(numMics,sf.numTimeFrames,sf.specLen);
sf.psi = cell(numMics,sf.numTimeFrames);
for pp = 1:numMics
    A0_stft{pp} = getForwardSTFT(ai{pp}(:,1), sf.window, sf.noverlap, sf.nfft);
    A1_stft{pp} = getForwardSTFT(ai{pp}(:,2), sf.window, sf.noverlap, sf.nfft);
    A2_stft{pp} = getForwardSTFT(ai{pp}(:,3), sf.window, sf.noverlap, sf.nfft);
    A3_stft{pp} = getForwardSTFT(ai{pp}(:,4), sf.window, sf.noverlap, sf.nfft);
    for ii = 1:sf.numTimeFrames
        A_temp = [A0_stft{pp}(:,ii) A1_stft{pp}(:,ii) A2_stft{pp}(:,ii) A3_stft{pp}(:,ii)];
        [r_I{pp,ii}, sf.psi{pp,ii}] = merimaa2005(A_temp);
        for kk = 1:sf.specLen
            % DOA for mic p at each time-frequency bin
            s_p{pp,ii,kk} = r_I{pp,ii}(kk,:) / norm(r_I{pp,ii}(kk,:));
        end
    end
    sf.p{pp} = A0_stft{pp} * sqrt(4*pi); % Pressure at Each Mic
end

%% Estimate Source Positions
sf.s_0 = cell(sf.numTimeFrames,sf.specLen);
for ii = 1:sf.numTimeFrames
    for kk = 1:sf.specLen
        % Triangulated source position for each time-frequency bin
        sf.s_0{ii,kk} = AmbiNav_TriangulateSource(u, s_p(:,ii,kk), force2D);
    end
end

end

function [rI, D] = merimaa2005(B)

rho_0 = 1.225; % density of air in kg/m^3
Z0 = rho_0 * getSoundSpeed(); % acoustic impedance of air in kg/(m^2 s)

W = B(:,1)/sqrt(2);
Y = B(:,2)/sqrt(3);
Z = B(:,3)/sqrt(3);
X = B(:,4)/sqrt(3);

ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

Xp = X*ex + Y*ey + Z*ez;

temp = real(diag(conj(W))*Xp);
rI = (sqrt(2)/Z0)*temp;

F = sqrt(dot(temp,temp,2));
E = abs(W).^2 + dot(Xp,Xp,2)/2;
D = 1 - sqrt(2)*F./E;

end
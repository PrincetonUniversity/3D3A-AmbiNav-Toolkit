function a_out = thiergart2013_synthesis(sf, u, Lo, r, refMic)

%% Estimate Direct and Diffuse Pressures at IPLS
S_dir  = zeros(sf.specLen,sf.numTimeFrames);
S_diff = zeros(sf.specLen,sf.numTimeFrames);
for ii = 1:sf.numTimeFrames
    for kk = 1:sf.specLen
        d_source = sf.s_0{ii,kk} - u{refMic};
        S_ref = sf.p{refMic}(kk,ii); % Pressure at Reference Mic
        
        Gamma_ref = (1/sf.psi{refMic,ii}(kk)) - 1; % Eqs. (26) and (27)
        G_dir_ref = sqrt(Gamma_ref / (1 + Gamma_ref)); % Eq. (8)
        S_dir_ref = G_dir_ref * S_ref; % Eq. (7)
        H_dir = conj(ambPointSource(0,sf.kVec(kk),norm(d_source))); % = conj((1/norm(d_source)) * exp(1i*sf.kVec(kk)*norm(d_source)));
        S_dir(kk,ii) = S_dir_ref ./ H_dir; % Eq. (6)
        
        G_diff_ref = sqrt(1 / (1 + Gamma_ref)); % Eq. (14a)
        S_diff_ref = G_diff_ref * S_ref; % Eq. (13)
        S_diff(kk,ii) = S_diff_ref ./ 1; % Eq. (12), with H_diff = 1
    end
end

%% Assemble Output Signals
outLen = (sf.numTimeFrames - 1)*(sf.nfft - sf.noverlap) + sf.nfft;
No = (Lo + 1)^2;
a_out = zeros(outLen,No);
for nn = 1:No
    [l,~] = getAmbOrder(nn-1);
    Q_factor = sqrt(1 / (4*pi));
    B_n = zeros(sf.nfft,sf.numTimeFrames);
    for ii = 1:sf.numTimeFrames
        for kk = 1:sf.specLen
            d_source = sf.s_0{ii,kk} - r;
            Y = AmbiNav_SphericalHarmonicY(l, d_source);
            C_factor = Y(nn);
            H_dir = conj(ambPointSource(l,sf.kVec(kk),norm(d_source)));
            B_n(kk,ii) = (C_factor * H_dir * S_dir(kk,ii)) + (Q_factor * 1 * S_diff(kk,ii)); % Eq. (33), with H_diff = 1
        end
    end
    a_out(:,nn) = AmbiNav_InverseSTFT(B_n, sf.window, sf.noverlap, sf.nfft);
end

end

function d = ambPointSource(l,k,r)

if k(1)==0
    DIM = find(size(k) ~= 1, 1, 'first');
    k = k(2:end);
    dropZero = true;
    zeroVal = +~l;
else
    dropZero = false;
end

d = (1i^(l+1)*k).*AmbiNav_SphericalHankelH(l,1,k*r);

if dropZero
    if ~isempty(DIM)
        d = cat(DIM,zeroVal,d);
    elseif dropZero && isempty(DIM)
        d = zeroVal;
    end
end

end
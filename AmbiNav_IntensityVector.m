function rI = AmbiNav_IntensityVector(B)
%MERIMAA2005_I Merimaa's intensity vector from ambisonics signals.
%   I = MERIMAA2005_I(X) computes the intensity vector using ambisonics
%   signals X. The signals should be in columns in ACN order. I will be a
%   K-by-3 matrix, where K is the number of rows in X.
%   
%   Following Merimaa and Pulkki (2005) Spatial Impulse Response Rendering
%       I: Analysis and Synthesis.

narginchk(1, 1);

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

rI = (sqrt(2)/Z0)*real(diag(conj(W))*Xp);

end
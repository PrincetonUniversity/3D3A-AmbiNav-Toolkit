function D = AmbiNav_DiffusenessParameter(B)
%MERIMAA2005_D Merimaa's diffuseness parameter from ambisonics signals.
%   D = MERIMAA2005_D(X) computes the diffuseness parameter using
%   ambisonics signals X. The signals should be in columns in ACN order. D
%   will be a K-by-1 vector, where K is the number of rows in X.
%   
%   Following Merimaa and Pulkki (2005) Spatial Impulse Response Rendering
%       I: Analysis and Synthesis.

narginchk(1, 1);

W = B(:,1)/sqrt(2);
Y = B(:,2)/sqrt(3);
Z = B(:,3)/sqrt(3);
X = B(:,4)/sqrt(3);

ex = [1 0 0];
ey = [0 1 0];
ez = [0 0 1];

Xp = X*ex + Y*ey + Z*ez;

temp = real(diag(conj(W))*Xp);
F = sqrt(dot(temp,temp,2));
E = abs(W).^2 + dot(Xp,Xp,2)/2;
D = 1 - sqrt(2)*F./E;

end
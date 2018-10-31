function val = AmbiNav_SphericalHankelH(n,k,x)

norm = sqrt(pi./(2*x));
sgn = 2*(x>=0)-1;
%val = norm.*sgn.*besselh(n+0.5,k,x);
val = norm.*sgn.*(besselj(n+0.5,x)-((-1)^k).*1i.*bessely(n+0.5,x));

end
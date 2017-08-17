function U=Green_Grad(k,z,z0)
%% The calculation of the gradient of Green function.
z1=z(1,:)-z0(1,:);
z2=z(2,:)-z0(2,:);
d=sqrt(z1.*z1+z2.*z2);
val  = besselh(1,k*d);
if d==0
    U=[0;0];
else
    U = 1i/4*k*[val.*z1./d; val.*z2./d];
end


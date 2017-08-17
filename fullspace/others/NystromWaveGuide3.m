function U=NystromWaveGuide3(n, n_src, n_recv, bctype,wavenumber, source,receiver, H)
%%  source: positon of source
%% receive: position of receive

%% compute the quadrature weights;

w = quad_weights(n);
R = zeros(2*n);

for k=1:2*n
     idx=[k:2*n];
     R(idx,k)=w([1:2*n-k+1]);
     R(k,k)=R(k,k)/2;  %% for convinience
end
R=(R+R');
%% discrete point
node = 0:2*n-1;
t = pi*node(:)/n;
if bctype==1
    [x1,x2]=circlebc(t,1);   
    [dx1,dx2]=circlebc(t,2);
    [ddx1,ddx2]=circlebc(t,3);
else
    [x1,x2]=kite(t,1);   
    [dx1,dx2]=kite(t,2);
    [ddx1,ddx2]=kite(t,3);
end

r = zeros(2*n);

for k=1:2*n
    for j=1:2*n
        r(k,j)=sqrt((x1(j)-x1(k))* (x1(j)-x1(k))+(x2(j)-x2(k))* (x2(j)-x2(k)));
    end
end

L1 = zeros(2*n);
M1 = zeros(2*n);
L2 = zeros(2*n);
M2 = zeros(2*n);
L = zeros(2*n);
Ceuler = 0.577215665;
distance = sqrt( dx1.*dx1+dx2.*dx2 );

for j=1:2*n
    for k=1:2*n
        if (j==k)
            dist = distance(k);
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            v1=[x1(j);x2(j)]; 
            v2=[x1(k);x2(k)];
            M2(j,k) = (1i/2-Ceuler/pi-1/(2*pi)*log(wavenumber^2/4*dist*dist))*dist;
            L(j,k) = ( ( -1/pi*(log(pi/(4*H))-log(tan(pi/H*x2(j)/2))) - (1i/2-Ceuler/pi-1/(pi)*log(wavenumber/2))) + 2*Green2DWaveGuide(wavenumber,H,400,6,v1,v2) )*dist;
        else
            dist = distance(k);
           
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            M2(j,k) = 1i/2*besselh(0,wavenumber*r(j,k))*dist - M1(j,k)*log(4*sin((t(j)-t(k))/2)^2);
            
            v1=[x1(j);x2(j)]; 
            v2=[x1(k);x2(k)];
            L(j,k) =  2*Green2DWaveGuide(wavenumber,H,400,5,v1,v2)*distance(k);
        end
       
    end
end
%% the linear System is A = I-(L1+ik*M1 +pi/n*(L2+ik*M2))



%% the right hand side is double of incidient 

%% 
for is=1:n_src
    for ix=1:2*n
        v1 = [x1(ix);x2(ix)];
        f(ix,is) = -2*Green2DWaveGuide( wavenumber, H, 400, 4, v1, source(:,is) );
    end
end

 S = (R.*M1+pi/n*M2) + pi/n*L;
%% the solution is the density function on the boundary of D
 %% phi = A\f;
 phi = S\f;

%% Composite Trapzitol Formula for Computing Near Field scattering data u_s




U = zeros(n_recv,n_src);

for k=1:n_recv
    for ix=1:2*n
         v1 = [x1(ix);x2(ix)];
         g(ix) = Green2DWaveGuide( wavenumber, H, 400, 4, v1, receiver(:,k) );
    end
    for j=1:n_src
        temp =  (  (distance').*g )*phi(:,j);
        U(k,j)=pi/n*temp;
    end
end
function U = NystromTransmission(n, n_src, n_recv, bctype,wavenumber, wavenumber2, source,receiver)
%% Nystrom methods for Transmission problem in waveguide setting.


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

%% quadrature weirghts for T operator;


%% Nystrom Methods : Assembling Integral Kernel Matrix %%

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




S1 = layerpoteintial(wavenumber, n, bctype, 1);
S2 = layerpoteintial(wavenumber2, n, bctype, 1);

K1 = - eye(2*n) +  layerpoteintial(wavenumber, n, bctype, 2);
K2 =  eye(2*n) + layerpoteintial(wavenumber2, n, bctype, 2);

Mat=[S1 -S2; K1 -K2];

distance = sqrt( dx1.*dx1+dx2.*dx2 );


%% the right hand side is double of incidient 
%% outer unit normal direction.
n_x = dx2(:)./distance(:);
n_z = -dx1(:)./distance(:);

f = zeros(4*n,n_src);
for is=1:n_src
    
    g = Green(wavenumber,[x1';x2'],source(:,is)*ones(1,2*n));
    gd = Green_Grad(wavenumber,source(:,is)*ones(1,2*n),[x1';x2']);
    f(1:2*n,is) = - 2 * g(:)  ;
    f(2*n+1:4*n,is) = - 2 * ( n_x(:).*( gd(1,:).') + n_z(:).* (gd(2,:).') );
end


%% Solve the linear system for the density function in the outer domain.
phi=Mat\f;
%% the solution is the potential on the boundary of D
%% phi = A\f;

%% Composite Trapzitol Formula for Computing Far fields pattern
U = zeros(n_recv,n_src);

for k=1:n_recv
   g = Green(wavenumber,receiver(:,k)*ones(1,2*n),[x1';x2']);
    for j=1:n_src
        temp =  (  (distance').*g )*phi(1:2*n,j);
     %%   temp = (distance'.*g1)*phi(:,j);
        U(k,j)=pi/n*temp;
    end
end

end
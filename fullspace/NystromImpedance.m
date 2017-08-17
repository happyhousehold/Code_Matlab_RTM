function U = NystromImpedance(n, n_src, n_recv, bctype, wavenumber, source, receiver)
%% Nystrom methods for Impedent problem
%% bctype==0 for Dirichlet BC; bctype==1 for Impenetrable BC


%% U_s  = \int_{\gamma_{D}} G(x,y)\phi(y) ds_{y}) 
%% the jump relation dU/dn = (-1/2 + K'+i k lambda S)\phi = -(du_{i}/dn+iklambda*u_{i})


%% Incidient direction: d=(cos(theta),sin(theta))


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

%% quadrature weirghts for T operator;


%% Nystrom Methods : Assembling Integral Kernel Matrix 

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



M1 = zeros(2*n);
M2 = zeros(2*n);
H1 = zeros(2*n);
H2 = zeros(2*n);

Ceuler = 0.577215665;
distance = sqrt( dx1.*dx1+dx2.*dx2 );

for j=1:2*n
    for k=1:2*n
        if j==k
            dist = distance(k);
            H1(j,k) = 0;
            H2(j,k) = 1/(2*pi)*( - dx1(j)*ddx2(j) + dx2(j)*ddx1(j))/dist;
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            M2(j,k) = (i/2-Ceuler/pi-1/(2*pi)*log(wavenumber^2/4*dist*dist))*dist;
            

        else
            dist = distance(k);
            lg4s = log(4*sin((t(j)-t(k))/2)^2);
            temp = dx2(j)*(x1(k)-x1(j))-dx1(j)*(x2(k)-x2(j));
            dxdx = dx1(j)*dx1(k)+dx2(j)*dx2(k) ;
          
            H1(j,k) = -wavenumber/(2*pi)*temp*besselj(1,wavenumber*r(j,k))/r(j,k)*dist;
            H2(j,k) = i*wavenumber/2*temp*besselh(1,wavenumber*r(j,k))/r(j,k)*dist - H1(j,k)*lg4s;
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            M2(j,k) = i/2*besselh(0,wavenumber*r(j,k))*dist - M1(j,k)*log(4*sin((t(j)-t(k))/2)^2);

        end
    end
end

%% the linear System is A = I-(L1+ik*M1 +pi/n*(L2+ik*M2))
 lambda = zeros(2*n,1);
%  lambda = 2 + x1.^3 + x2.^3; 
 lambda(1:n) = 1;
 lambda(n+1:end) = 1;

%  lambda(1:2*n) = 10000;
%% A = T + R.*K1 + pi/n*K2 + i*eta*diag(distance);
S = R.*M1 + pi/n*M2;
K = R.*H1 + pi/n*H2;
A =  diag(distance) - ( K + 1i*wavenumber*diag(distance(:).*lambda(:))*S);


%% the right hand side is double of incidient 
for is=1:n_src
    g1  = Green(wavenumber, source(:,is)*ones(1,2*n),     [x1';x2']);
    gd1 = Green_Grad(wavenumber,source(:,is)*ones(1,2*n), [x1';x2']);
    f(:,is) = 2*(  dx2(:).*( gd1(1,:).') - dx1(:).* (gd1(2,:).')  + 1i*wavenumber*( distance(:).*lambda(:).*g1(:) ) );
end


%% the solution is the potential on the boundary of D
 phi = A\f;
%phi = S\f;
%% Composite Trapzitol Formula for Computing Far fields pattern

for j=1:n_src
    for k=1:n_recv
       g1 = Green(wavenumber,receiver(:,k)*ones(1,2*n),[x1';x2']);
        temp =  (distance'.*g1)*phi(:,j);
       %% temp = (wavenumber*( xhat*dx2 - yhat*dx1) + eta*distance).*(exp(-i*wavenumber*(xhat*x1+yhat*x2))).*phi(:,j);
       U(k,j)=pi/n*sum(temp);
    end
end

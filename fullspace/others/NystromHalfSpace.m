function U=NystromHalfSpace(n, n_src, n_recv, bctype,wavenumber, source,receiver)
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
          
            M2(j,k) = (i/2-Ceuler/pi-1/(2*pi)*log(wavenumber^2/4*dist*dist))*dist;
        else
            dist = distance(k);
            temp = dx2(k)*(x1(j)-x1(k))-dx1(k)*(x2(j)-x2(k));
           
            M1(j,k) = -1/(2*pi)*besselj(0,wavenumber*r(j,k))*dist;
            M2(j,k) = i/2*besselh(0,wavenumber*r(j,k))*dist - M1(j,k)*log(4*sin((t(j)-t(k))/2)^2);
        end
        D = sqrt( (x1(j)-x1(k)).^2 + (x2(j)+x2(k)).^2 );
        L(j,k) = - i/2*besselh(0,wavenumber*D)*distance(k);
    end
end
%% the linear System is A = I-(L1+ik*M1 +pi/n*(L2+ik*M2))




%% the right hand side is double of incidient 

%% 
for is=1:n_src
    f(:,is) = -2*(Green(wavenumber, source(:,is)*ones(1,2*n),[x1';x2']) - Green(wavenumber, source(:,is)*ones(1,2*n),[x1';-x2']));
end

 S = (R.*M1+pi/n*M2) + pi/n*L;
%% the solution is the potential on the boundary of D
 %% phi = A\f;
 phi = S\f;

%% Composite Trapzitol Formula for Computing Far fields pattern




U = zeros(n_recv,n_src);
for j=1:n_src
    for k=1:n_recv
       
        g1 = Green(wavenumber,receiver(:,k)*ones(1,2*n),[x1';x2']);
        g2 = Green(wavenumber,receiver(:,k)*ones(1,2*n),[x1';-x2']);
        temp =  (  (distance').*(g1-g2) )*phi(:,j);
     %%   temp = (distance'.*g1)*phi(:,j);
        U(k,j)=pi/n*temp;
    end
end

% Nx = 201;
% Nz = 201;
% Irtm = zeros(Nx,Nz);
% x = linspace( -8, 8, Nx);
% z = linspace( -8, 8, Nz);
% Nfreq = 1;size(U)
% for ix = 1 : Nx
%     for iz = 1: Nz
%         for k=1:Nfreq
%             gr = conj(((Green(wavenumber,receiver,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)]))));
%             gs = conj(((Green(wavenumber,source,[x(ix)*ones(1,n_src); z(iz)*ones(1,n_src)]))));
%             Irtm(ix,iz)=Irtm(ix,iz)+sum( gs.* ( gr*U(:,:) ) );
%         end
%     end
% end
% 
% 
% figure
% [xx,zz]=meshgrid(x,z);
% subplot(1,3,1)
% imagesc(abs(imag(Irtm))');
% subplot(1,3,2);
% mesh(xx,zz,abs(imag(Irtm))');
% subplot(1,3,3)
% plot(x,imag(Irtm(101,:)));
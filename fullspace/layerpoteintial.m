function S=layerpoteintial(freq, n, bctype, type)
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


% [x1,x2]=kite(t,1);   
% [dx1,dx2]=kite(t,2);
% [ddx1,ddx2]=kite(t,3);

[x1,x2]=circlebc(t,1);   
[dx1,dx2]=circlebc(t,2);
[ddx1,ddx2]=circlebc(t,3);
% [dddx1,dddx2]=circlebc(t,4);

r = zeros(2*n);

for k=1:2*n
    for j=1:2*n
        r(k,j)=sqrt((x1(j)-x1(k))* (x1(j)-x1(k))+(x2(j)-x2(k))* (x2(j)-x2(k)));
    end
end


M1 = zeros(2*n);

M2 = zeros(2*n);



Ceuler = 0.577215665;
distance = sqrt( dx1.*dx1+dx2.*dx2 );
if type==1
    for j=1:2*n
        for k=1:2*n
            if (j==k)
                dist = distance(k);
                
                M1(j,k) = -1/(4*pi)*dist;
                M2(j,k) = 1/2*( i/2 - Ceuler/pi - 1/pi*log(freq*dist/2) )*dist;
                
            else
                dist = distance(k);
                lg4s = log(4*sin((t(j)-t(k))/2)^2);
                
                M1(j,k) = -1/(4*pi)*besselj(0,freq*r(j,k))*dist;
                M2(j,k) = i/4*besselh(0,freq*r(j,k))*dist - M1(j,k)*lg4s;
                
            end
        end
    end
elseif type==2
    for j=1:2*n
        for k=1:2*n
            if (j==k)
                dist = distance(k);
                M1(j,k) = 0;
                M2(j,k) = 1/(4*pi)*( dx2(j)*ddx1(j) - dx1(j)*ddx2(j) )/dist;
                
            else
                dist = distance(k);
                lg4s = log(4*sin((t(j)-t(k))/2)^2);
                temp = distance(k)/distance(j)*( dx2(j)*(x1(j)-x1(k))-dx1(j)*(x2(j)-x2(k)) );
                
                M1(j,k) = freq/(4*pi)*temp*besselj(1,freq*r(j,k))/r(j,k);
                M2(j,k) = -i*freq/4*temp*besselh(1,freq*r(j,k))/r(j,k) - M1(j,k)*lg4s;
                
            end
        end
    end
end

S = R.*M1 + pi/n*M2;


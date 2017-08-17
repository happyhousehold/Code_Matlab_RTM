function w=quad_weights(n)
%% w is of length 2n-1
w = zeros(2*n,1);
pts = [1:n-1];
for k=0:2*n-1
    temp = cos(k*pi*pts/n)./pts;
%%    w(k+1) = -2*pi/n*sum(temp)-(-1)^k*pi/(n*n); 
    w(k+1) = -2*pi/n*(sum(temp)+(-1)^k/(2*n));
end
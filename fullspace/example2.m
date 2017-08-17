%% ---------------- Parameters Setting---------------------------------------------%%
%clear;
bctype = 1;  %% Types of obstacle: bctype=1 for Circle; bctype=2 for Kite;  bctype=3 for Rounded-Square  

n_src = 256; %% number of source and receiver
n_recv = 256;

npts = 64 ; %% discretization point of integral equation

freq = 2*pi*2; 


source   = zeros(2,n_src);
receiver = zeros(2,n_recv);

theta_r = (0:n_recv-1)*2*pi/n_recv;
theta_s = (0:n_src-1)*2*pi/n_src ;
R=15;
source(1,:)   = R*cos(theta_s); source(2,:)   = R*sin(theta_s);
receiver(1,:) = R*cos(theta_r); receiver(2,:) = R*sin(theta_r);

%% ---------------- End of Parameters Setting----------------------------------------%%

%% ----------------Nystrom Method for Synthesizing Scattering Data ------------------%%
tic
U = Nystrom_MultipleObjects(npts, n_src, n_recv, 1,freq, source,receiver); %% multiple obstacles
%% U = Nystrom_MultipleObjects2(npts, n_src, n_recv, 1,freq, source,receiver); %% multiple obstacles
disp('Nystrom Method for Synthesizing  Obstacle Scattering Data Complete......');

toc

%% ----------------End of  Synthesizing Scattering Data -----------------------------%%

%% --------------------Reverse Time Migration for Obstacle Imaging-------------------%%
%% Sampling domain;
Nx = 201;
Nz = 201;
x = linspace( -10, 10, Nx);
z = linspace( -10, 10, Nz);

flag=1;
Irtm =  RTMFullAperture(n_src, n_recv, U, freq, source, receiver, x, z, R, flag);

fprintf('The maximum value of imag(Irtm) is  %f\n',max(max(imag(Irtm))));
fprintf('The minimum value of imag(Irtm) is  %f\n',min(min(imag(Irtm))));

figure, 
subplot(1,2,1);imagesc(x,z,real(Irtm)');colorbar; set(gca, 'YDir', 'normal');
subplot(1,2,2);imagesc(x,z,imag(Irtm)');colorbar; set(gca, 'YDir', 'normal');
figure,
[xx,zz]=meshgrid(x,z); 
subplot(1,2,1);mesh(xx,zz,(real(Irtm))');
subplot(1,2,2);mesh(xx,zz,(imag(Irtm))');
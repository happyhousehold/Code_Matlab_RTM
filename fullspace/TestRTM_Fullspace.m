%% ---------------- Parameters Setting---------------------------------------------%%
clear;
bctype = 1;  %% Types of obstacle: bctype=1 for Circle; bctype=2 for Kite;  bctype=3 for Rounded-Square  

n_src = 64; %% number of source and receiver
n_recv = 64;

npts = 64 ; %% discretization point of integral equation

numd=[2];
freq = 2*pi*numd; 
nf=length(numd);

flag=1;
source   = zeros(2,n_src);
receiver = zeros(2,n_recv);

theta_r = (0:n_recv-1)*2*pi/n_recv;
theta_s = (0:n_src-1)*2*pi/n_src ;
R=10;
source(1,:)   = R*cos(theta_s); source(2,:)   = R*sin(theta_s);
receiver(1,:) = (R+3)*cos(theta_r); receiver(2,:) = (R+3)*sin(theta_r);

%% ---------------- End of Parameters Setting----------------------------------------%%

%% ----------------Nystrom Method for Synthesizing Scattering Data ------------------%%
tic
for k=1:nf
    if flag==1
          U(:,:,k)=Nystrom_NearFieds(npts, n_src, n_recv, bctype,freq(k),source,receiver);   %% near fields pathern
        %% U(:,:,k) = NystromImpedance(npts, n_src, n_recv, bctype,freq(k), source, receiver); % Nystrom methods for Impedent problem
        %% U(:,:,k) = NystromTransmission(npts, n_src, n_recv, bctype,freq(k), freq(k)/2, source,receiver); %% transmission data
         %% U(:,:,k) = Nystrom_MultipleObjects2_pl(npts, n_src, n_recv, 1,freq(k), source,receiver); %% multiple obstacles
    else
            U(:,:,k)=Nystrom(npts, n_src, n_recv, bctype,freq(k));   %% far fields pathern
    end
end
disp('Nystrom Method for Synthesizing  Obstacle Scattering Data Complete......');

toc

%% ----------------End of  Synthesizing Scattering Data -----------------------------%%

%% --------------------Reverse Time Migration for Obstacle Imaging-------------------%%
%% Sampling domain;
Nx = 201;
Nz = 201;

if bctype==1
    x = linspace( -3, 3, Nx);
    z = linspace( -3, 3, Nz);
elseif bctype==2
    x=linspace(-3,3,Nx);
    z=linspace(1.5,6.5,Nz);
end
 [xx,zz]=meshgrid(x,z); 

Irtm =  RTMFullAperture(n_src, n_recv, U, freq, source, receiver, x, z, R, flag);

if bctype==1
    disp('The obstacle we are imaging is a circle.');
elseif bctype==2
    disp('The obstacle we are imaging is a kite.');
end

fprintf('The maximum value of imag(Irtm) is  %f\n',max(max(imag(Irtm))));
fprintf('The minimum value of imag(Irtm) is  %f\n',min(min(imag(Irtm))));

figure,
imagesc(x,z,imag(Irtm)');colormap('jet');colorbar; set(gca, 'YDir', 'normal');
figure,
mesh(xx,zz,(imag(Irtm))');colormap('jet');colorbar; set(gca, 'YDir', 'normal');

% figure, 
% subplot(1,2,1);imagesc(x,z,real(Irtm)');colorbar; set(gca, 'YDir', 'normal');
% subplot(1,2,2);imagesc(x,z,imag(Irtm)');colorbar; set(gca, 'YDir', 'normal');
% figure,

% subplot(1,2,1);mesh(xx,zz,(real(Irtm))');
% subplot(1,2,2);mesh(xx,zz,(imag(Irtm))');
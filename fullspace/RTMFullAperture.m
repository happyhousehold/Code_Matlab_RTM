function Irtm =  RTMFullAperture(n_src, n_recv, U, freq, source, receiver, x, z, R, flag)
%% Reverse Time Migration for obstacle scattering data available on a circle domain
%% n_src, n_recv denotes the number of source and receiver
%% U(n_recv,n_src,freq)
Nfreq = length(freq);
Nx = length(x);
Nz = length(z);
Irtm = zeros(Nx,Nz);
tic

if flag==1
    %% Near fields RTM imaging
    scaling = (2*pi*R)^2/n_src/n_recv;
%     scaling = 1/Nfreq*(receiver(end,1)-receiver(1,1))*(source(end,1)-source(1,1))/n_src/n_recv;
    for ix=1:Nx
        for iz=1:Nz
            for k=1:Nfreq
                if n_src == n_recv
                    gs =  conj(Green(freq(k),source,[x(ix)*ones(1,n_src); z(iz)*ones(1,n_src)])) ;
                    gr = gs;
                else
                    gr = conj(Green(freq(k),receiver,[x(ix)*ones(1,n_recv); z(iz)*ones(1,n_recv)])) ;
                    gs = conj(Green(freq(k),source,[x(ix)*ones(1,n_src); z(iz)*ones(1,n_src)])) ;
                end
                    Irtm(ix,iz)=Irtm(ix,iz) + freq(k)^2*scaling*sum( gs.* ( gr*U(:,:,k) ) ); %% correlation imaging condition.
                    
                    
               %% Irtm(ix,iz)=Irtm(ix,iz) + scaling*freq(k)*sum(  ( gr*U(:,:,k) )./conj(gs) ); %% deconvolution imaging condition
             %%   Irtm(ix,iz)=Irtm(ix,iz) + scaling*freq(k)^2*sum( abs( ( gr*U(:,:,k) ) ) ); %% backpropagation imaging
         
            end
        end
    end
else
    scaling = 1/n_src/n_recv;
    for ix=1:Nx
        for iz=1:Nz
            for k=1:Nfreq
                coef = exp(1i*pi/4)/sqrt(8*pi*freq(k));
                gr = exp( 1i*freq(k)*( receiver(1,:)*x(ix)+receiver(2,:)*z(iz) ) );
                gs = exp( -1i*freq(k)*( source(1,:)*x(ix)+source(2,:)*z(iz) ) );
                Irtm(ix,iz)=Irtm(ix,iz) + 1/4*1/coef*scaling*sum( gs.* ( gr*U(:,:,k) ) );
            end
        end
    end
end
disp('Imaging via RTM method is completed......');
toc


    
            
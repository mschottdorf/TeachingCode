close all;
clear all;

rng(1);

% Key parameters
g = 0.98;
sigma = 1.7;

nHyper = 22;
N = 512;
pixels_per_hc = N/nHyper;

% For symmetry breaking + time steps.
r = 0.01;  % Characterisitc timescale is 1/r
dt=0.1;
t=0;

KK=zeros(N,1);
for k=1:N
  KK(k)=(2*(k-1)/N-1)*N/(2*nHyper);
end
reg_sh_kernel=zeros(N);
for k=1:N
  for l=1:N
    reg_sh_kernel(k,l)=1./(1-dt*(r-(1-(KK(k)^2+KK(l).^2)).^2));
  end
end


ss = sigma*(N/nHyper);
h = zeros(N,N);
for x=1:N
    for y=1:N
        h(x,y) = exp(-((x-N/2)^2+(y-N/2)^2)/(2*ss^2));
    end
end
gaussian = fftshift(h) ./ sum(h(:)); % Kernel


%initial conditions
z_map = 1e-2*randn(N) + 1i*1e-2*randn(N);

tic
while(t<1e6)  
    
    do_plot = mod(round(t),100)==0;
    if do_plot 
        figure(2),clf;
        subplot(1,3,1)
        imagesc(angle(z_map));
        colormap(hsv)
        freezeColors;
        title("Current OPM")

        subplot(1,3,2)
        imagesc(KK,KK,abs(fftshift(fft2(z_map-mean(mean(z_map))))).^2);
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        colormap(flipud(gray))
        freezeColors;
        title("Current OPM spectrum")

        subplot(1,3,3)
        hist(abs(z_map(:)))
        xlim([0, 1])
        title("Amplitudes")
        drawnow;
    end


    N1 = ifft2( fft2(gaussian) .* fft2(abs(z_map).^2) ) .* z_map;
    N2 = ifft2( fft2(gaussian) .* fft2(z_map.^2) ) .* conj(z_map)./2;

    nonlocal = dt*(g-2)*(N1+N2);
    local = dt*(1-g)*(abs(z_map).^2).*z_map;
    
    % Update z
    u=fftshift(fft2(z_map));
    z_map = ifft2(ifftshift(u.*reg_sh_kernel)) + local + nonlocal;
    t=t+dt;
    disp([t, mean(abs(z_map(:)))]);

end
toc



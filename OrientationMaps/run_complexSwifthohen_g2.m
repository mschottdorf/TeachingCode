close all;
clear all;

rng(1);

r = 0.01;
g = 2;

nHyper = 22;
N = 512;
pixels_per_hc = N/nHyper;

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

%initial conditions
z_map = 1e-2*randn(N) + 1i*1e-2*randn(N);

tic
while(t<1e4)  % 10h

    if mod(round(t),100)==0; 
        figure(22),clf;
        subplot(1,2,1)
        imagesc(angle(z_map));
        colormap(hsv)
        freezeColors;
        title("Current OPM")

        subplot(1,2,2)
        imagesc(KK,KK,abs(fftshift(fft2(z_map-mean(mean(z_map))))).^2);
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        colormap(flipud(gray))
        freezeColors;
        title("Current OPM spectrum")
  
        drawnow;
    end

    % Update z
    u=fftshift(fft2(z_map));

    z_map = ifft2(ifftshift(u.*reg_sh_kernel)) + dt*(1-g)*(abs(z_map).^2).*z_map;
    t=t+dt;
    disp([t, mean(abs(z_map(:)))]);
end
toc


figure(22),clf;
subplot(1,2,1)
imagesc(angle(z_map));
colormap(hsv)
freezeColors;
title("Current OPM, g=2 100000t")

subplot(1,2,2)
imagesc(KK,KK,abs(fftshift(fft2(z_map-mean(mean(z_map))))).^2);
xlim([-1.5,1.5])
ylim([-1.5,1.5])
colormap(flipud(gray))
freezeColors;
title("Current OPM spectrum, g=2 100000t")

drawnow;

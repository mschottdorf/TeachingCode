close all;
clear all;

rng(1);

r = 0.001;
kappa = 0.0001* (sqrt(r).^3) * 5;
gamma = (sqrt(r)/2) * (sqrt(3)/2.) * 0.5;

nHyper = 22;
N = 512;
pixels_per_hc = N/nHyper;

dt=0.1;
t=0;

%this is for the external field:
Nmodes = 18; % modes on the critical circle
l = randi([0 1],Nmodes,1)*2-1;
phase = rand(Nmodes,1) + 1i*rand(Nmodes,1);
z = create_ECP([N,N], nHyper, l, phase);
ext = abs(z).^2;


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

%%% Plotting external field
figure(1),clf;
subplot(1,2,1)
imagesc(angle(z));
colormap hsv;
freezeColors;
title("Bias field")

subplot(1,2,2)
imagesc(KK,KK,abs(fftshift(fft2(z-mean(mean(z))))).^2);
xlim([-1.5,1.5])
ylim([-1.5,1.5])
colormap(flipud(gray))
title("Bias spectrum")


%initial conditions
co_map = 1e-2*randn(N);

do_plot = false;
tic
while(t<1e5)  

    if do_plot 
        figure(2),clf;
        subplot(1,2,1)
        imagesc(co_map);
        colormap(gray)
        freezeColors;
        title("Current CO map")
        subplot(1,2,2)
        imagesc(KK,KK,abs(fftshift(fft2(co_map-mean(mean(co_map))))).^2);
        xlim([-1.5,1.5])
        ylim([-1.5,1.5])
        colormap(flipud(gray))
        freezeColors;
        title("Current CO spectrum")
    
        figure(3),clf;
        x=(-1.5:.01:1.5);
        plco = co_map-mean(mean(co_map));
        plco = plco./std(plco(:));
        z_ROI = ones(size(z));
        c_ROI = ones(size(plco));
        max_shift = (-min(x) + max(x))/2;
        shift_y = ((0:size(z,1)-1) - floor(size(z,1)/2))/pixels_per_hc;
        shift_x = ((0:size(z,2)-1) - floor(size(z,2)/2))/pixels_per_hc;
        num_pixels = xcorr2_fft(double(c_ROI),double(z_ROI));
        value = xcorr2_fft(co_map, abs(z).^2)./num_pixels;
        E = value(abs(shift_y)<=max_shift, abs(shift_x)<=max_shift);
        imagesc(x, x, E);
        colorbar;
        title(['Current alignment at t= ', num2str(t)])
    
        drawnow;
    end

    % Update CO
    u=fftshift(fft2(co_map));
    co_map = real(ifft2(ifftshift(u.*reg_sh_kernel)))  + dt*(gamma*co_map.^2 - co_map.^3 - kappa*ext);
    t=t+dt;
    disp(t)
end
toc


function data=create_ECP(dim,nHyper,l,phase)
    x=1:dim(1);
    y=1:dim(2);
    [Y,X]=meshgrid(y,x);
    phi=(0:length(l)-1)*pi/length(l);
    data=zeros(dim);
    k=2*pi*nHyper/dim(1);
    
    for i=1:length(l)
        if l(i)~=0
            data=data+exp(1i*k*l(i)*(cos(phi(i))*X+sin(phi(i))*Y)+1i*phase(i));
        end
    end
end
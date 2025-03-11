

beta = 1; % width of the spectrum. Larger -> narrower

z = create_OPM(beta, [512,512], 10);

[count,PWxList,PWyList,signList] = pw_finder_withsign(z);

figure(1)
imagesc(angle(z));
daspect([1 1 1])
colormap hsv;
hold on;
scatter(PWxList, PWyList, 'ko', 'filled')
title(count/(10*10))  % Notice convergence to pi for beta large

function psi=create_OPM(beta,dim,nHyper)

    kx=(-dim(1)/2:dim(1)/2-1)/nHyper;
    ky=(-dim(2)/2:dim(2)/2-1)/nHyper;
    [KX,KY]=meshgrid(ky,kx);
    
    Aln=log(2)+(1+beta)*gammaln((2+beta)/2)-(2+beta)*gammaln((1+beta)/2);
    Bln=2*gammaln((2+beta)/2)-2*gammaln((1+beta)/2);
    
    P_k=exp(Aln)*sqrt(KX.^2+KY.^2).^beta.*exp(-(KX.^2+KY.^2)*exp(Bln))/nHyper^2/2/pi;
    
    psi_fft=dim(1)*dim(2)*sqrt(P_k).*(randn(dim)+1i*randn(dim))/sqrt(2);
    psi=ifft2(fftshift(psi_fft));
end
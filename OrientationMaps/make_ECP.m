
N = 18; % modes on the critical circle

l = randi([0 1],N,1)*2-1;
phase = rand(N,1) + 1i*rand(N,1);

z = create_ECP([512,512], 22, l, phase);

[count,PWxList,PWyList,signList] = pw_finder_withsign(z);

figure(1)
imagesc(angle(z));
colormap hsv;
hold on;
scatter(PWxList, PWyList, 'ko', 'filled')
title(count/(22*22))


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
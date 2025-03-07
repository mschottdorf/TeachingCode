function C = xcorr2_fft(A,B)
    % Calculate 2D correlation of A and B in Fourier Domain, return the same
    % size as A.
    % Based on xcorr2_fft by Alessandro Masullo

    if nargin == 1
        B = A;
    end
    assert(isreal(A) && isreal(B),'Error in calculating cross correlation')

    % input dimensions
    adim = size(A);
    bdim = size(B);

    % Cross-correlation dimension
    cdim = adim+bdim-1;

    bpad = zeros(cdim);
    apad = zeros(cdim);

    apad(1:adim(1),1:adim(2)) = A;
    bpad(1:bdim(1),1:bdim(2)) = B(end:-1:1,end:-1:1);

    ffta = fft2(apad);
    fftb = fft2(bpad);
    C = real(ifft2(ffta.*fftb));

    % Return same dimensions as A
    C = C(adim(1)-floor(bdim(1)/2):adim(1)+ceil(bdim(1)/2-1),...
        adim(2)-floor(bdim(2)/2):adim(2)+ceil(bdim(2)/2)-1);
end
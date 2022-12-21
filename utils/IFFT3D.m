function res = IFFT3D(data)
    % 3D Fourier transform along the coil direction
    data = permute(data, [1 2 4 3]);
    for icoil = 1:size(data, 4)
        res(:,:,:,icoil) = ifft3c(data(:,:,:,icoil));
    end
    res = permute(res, [1 2 4 3]);
end

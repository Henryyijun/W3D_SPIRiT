
function res = FFT3D(data)
    data = permute(data, [1 2 4 3]);
    for icoil = 1:size(data, 4)
        res(:,:,:,icoil) = fft3c(data(:,:,:,icoil));
    end
    res = permute(res, [1 2 4 3]);
end
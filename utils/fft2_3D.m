function [ fcoe ] = fft2_3D( im )
%   fft2_3D         multiple 2D fft
%   input:
%   im              3D data 
%   output; 
%   fcoe            3D fourier coefficient
%
%   Channing Lau
%   Dec. 2017

fcoe = zeros(size(im));
for i = 1:size(im,3)
    %fcoe(:,:,i) = fftshift(fft2(im(:,:,i)));
    fcoe(:,:,i) = (fft2(im(:,:,i)));
end


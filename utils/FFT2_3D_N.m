function [ fcoe ] = FFT2_3D_N( Img )
%   FFT2_3D         multiple 2D fft with normalizaiton 
%   input:
%   Img              3D data 
%   output; 
%   fcoe            3D fourier coefficient
%
%  
%   Dec. 2019
global gpu_enable
if gpu_enable == true
	fcoe = zeros(size(Img), 'gpuArray');
else
  fcoe = zeros(size(Img));
end

for i = 1:size(Img,3)
    %fcoe(:,:,i) = fftshift(fft2(ifftshift(Img(:,:,i))));
    fcoe(:,:,i) = fftshift(fft2((Img(:,:,i))));
end

fcoe = fcoe/sqrt(size(Img,1)*size(Img,2));
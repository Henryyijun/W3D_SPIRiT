function [ Img ] = IFFT2_3D_N( Fcoe )
%   IFFT2_3D         multiple 2D ifft with normalizaiton 
%   input:
%    Fcoe            multiple 2D fourier coefficient
%   output; 
%
%    Img              3D data
%  
%   Dec. 2019
global gpu_enable
if gpu_enable == true
Img = zeros(size(Fcoe), 'gpuArray');
else
  Img = zeros(size(Fcoe));
end
  
for i = 1:size(Fcoe,3)
    %Img(:,:,i) = fftshift(ifft2(ifftshift(Fcoe(:,:,i))));
    Img(:,:,i) = (ifft2(ifftshift(Fcoe(:,:,i))));
end


Img = Img*sqrt(size(Fcoe,1)*size(Fcoe,2));
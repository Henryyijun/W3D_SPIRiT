function [ im ] = ifft2_3D( fcoe )
%   ifft2_3D         多张图片的二维反傅里叶变换
%   fcoe            傅里叶系数
%   im              三维数据
%
%   Channing Lau
%   Dec. 2017

im = zeros(size(fcoe));
for i = 1:size(fcoe,3)
    im(:,:,i) = ifft2(ifftshift(fcoe(:,:,i)));
end


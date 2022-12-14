function [ im ] = ifft2_3D( fcoe )
%   ifft2_3D         ����ͼƬ�Ķ�ά������Ҷ�任
%   fcoe            ����Ҷϵ��
%   im              ��ά����
%
%   Channing Lau
%   Dec. 2017

% im = zeros(size(fcoe));
for i = 1:size(fcoe,3)
    im(:,:,i) = ifft2(ifftshift(fcoe(:,:,i)));
end

im = im*sqrt(size(fcoe,1)*size(fcoe,2));
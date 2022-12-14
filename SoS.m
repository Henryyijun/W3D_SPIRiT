function [ im_SoS ] = SoS( im )
%   SoS         SoS(square root of squares) of a 3D image
%   input:
%   im          3D image
%   output:
%   im_SoS      SoS image
%
%   Channing Lau
%   Dec. 2017

im_SoS = sum(abs(im).^2,3).^0.5;

end


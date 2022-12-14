function [ mask,Real_R ] = FindMask( Size,R,AcsLoc )
%   FindMask            find the mask for the 3D data
%   input:              
%   Size                size of a image
%   R                   reduction factor 
%   AcsLoc              location of ACS lines
%   output:
%   mask                mask with 3D size
%   Real_R              total acceleration
%
%   Channing Lau
%   Dec. 2017

mask = zeros(Size);
if length(R) == 1
    SampleRow = mod(0:Size(1)-1,R);
    mask(SampleRow==0,:,:) = 1;
    mask(AcsLoc,:,:) = 1;
    Real_R = prod(Size)/sum(mask(:));
elseif length(R) == 2
    SampleRow = mod(0:Size(1)-1,R(1));
    SampleCol = mod(0:Size(2)-1,R(2));
    mask(SampleRow==0,SampleCol==0,:) = 1;
    mask(AcsLoc(:,1),AcsLoc(:,2),:) = 1;
    Real_R = prod(Size)/sum(mask(:));
else
	error('illegal input');
end
mask = logical(mask);

end


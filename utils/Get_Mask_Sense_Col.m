function mask_sense = Get_Mask_Sense(size_sense, row_size, col_size)
% Function: Get_Mask_Sense(size_sense, row_size, col_size) find the
%           lines that need to be filled, extend from the middle
%
% Parameter: size_sense: size of sense, 3D
%            row_size: size of row / 2-element vector
%            row_size: size of column
%
% Return: mask_sense: 0-1 matrix, show which k-space data are use to
%         estimate sensitivity
%
% Sunrise

mask_sense = zeros(size_sense(1), size_sense(2));
row_low = size_sense(1)/2 - row_size/2 + 1;
row_high = size_sense(1)/2 + row_size/2;
% 1. col_size is an 2-element vector, gc is the range of the ACS line(col)
if size(col_size, 2) == 2
    mask_sense(row_low:row_high, col_size(1):col_size(2)) = 1;
% 2. a random num of gc area
else
    row_low = size_sense(1)/2 - row_size/2 + 1;
    row_high = size_sense(1)/2 + row_size/2;
    col_low = size_sense(1)/2 - col_size/2 + 1;
    col_high = size_sense(1)/2 + col_size/2;
    mask_sense(row_low:row_high,col_low:col_high) = 1;
end

